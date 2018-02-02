/***************************************************************************
 *
 * Authors:    David Strelak (davidstrelak@gmail.com)
 *
 * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
 * 02111-1307  USA
 *
 *  All comments concerning this program package may be sent to the
 *  e-mail address 'xmipp@cnb.csic.es'
 ***************************************************************************/

#include "reconstruction_adapt_cuda/movie_alignment_correlation_gpu.h"


// FIXME: REMOVE
#include <sstream>
#include "data/filters.h"
#include "data/xmipp_fftw.h"
#define SSTR( x ) static_cast< std::ostringstream & >( \
        ( std::ostringstream() << std::dec << x ) ).str()
// FIXME

void ProgMovieAlignmentCorrelationGPU::loadFrame(const MetaData& movie, size_t objId, bool crop, Image<float>& out) {
	FileName fnFrame;
	movie.getValue(MDL_IMAGE, fnFrame, objId);
	if (crop) {
		Image<double>tmp;
		tmp.read(fnFrame);
		tmp().window(out(), yLTcorner, xLTcorner, yDRcorner, xDRcorner);
	} else {
		out.read(fnFrame);
	}
}


void ProgMovieAlignmentCorrelationGPU::loadData(const MetaData& movie,
		const Image<double>& dark, const Image<double>& gain,
		double targetOccupancy, const MultidimArray<double>& lpf) {
	Image<float> frame, gainF, darkF;
	MultidimArray<float> filter;
	// FIXME consider loading imgs in full size and cropping them on GPU
	bool cropInput = (yDRcorner != -1);
	// find input image size
	gainF.data.resize(gain());
	darkF.data.resize(dark());

	// FIXME extract
	Matrix1D<double> w(2);
	filter.initZeros(newYdim, newXdim/2+1);
	FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY2D(filter)
	{
		FFT_IDX2DIGFREQ(i, newYdim, YY(w));
		FFT_IDX2DIGFREQ(j, newXdim, XX(w));
		double wabs = w.module();
		if (wabs <= targetOccupancy)
			A2D_ELEM(filter,i,j) = lpf.interpolatedElement1D(
					wabs * newXdim);
	}


	loadFrame(movie, movie.firstObject(), cropInput, frame);
	int noOfImgs = nlast - nfirst + 1;
	size_t noOfFloats = noOfImgs * std::max(frame.data.yxdim, (frame.data.xdim/2+1) * frame.data.ydim * 2);
	float* imgs = new float[noOfFloats]();
	std::cout << "noOfFloats: " << noOfFloats << std::endl;
	int counter = -1;
	int paddedLineLength = (frame.data.xdim/2+1)*2;
	FOR_ALL_OBJECTS_IN_METADATA(movie) {
		counter++;
		if (counter < nfirst ) continue;
		if (counter > nlast) break;

		loadFrame(movie, __iter.objId, cropInput, frame);
		// FIXME optimize if necessary
		if (XSIZE(darkF()) > 0)
			frame() -= darkF();
		if (XSIZE(gainF()) > 0)
			frame() *= gainF();


//		************
//		IN-PLACE
//		************
		// copy line by line, adding offset at the end of each line
		// result is the same image, padded in the X dimension to (N/2+1)*2
		float* dest = imgs + ((counter-nfirst) * paddedLineLength * frame.data.ydim); // points to first float in the image
		for (int i = 0; i < frame.data.ydim; i++) {
			memcpy(dest + (paddedLineLength * i),
					frame.data.data + i*frame.data.xdim,
					frame.data.xdim * sizeof(float));
		}

//		************
//		OUT-OF-PLACE
//		************
//		// add image at the end of the stack (that is already long enough)
//		memcpy(imgs + ((counter-nfirst) * ((frame.data.xdim/2+1) * frame.data.ydim * 2)),
//				frame.data.data,
//				frame.data.yxdim * sizeof(float));
	}

	Image<float> aaaa((frame.data.xdim/2+1)*2, frame.data.ydim, 1, noOfImgs);
	aaaa.data.data = imgs;
	aaaa.write("images.vol");


//	float* result;
	size_t newFFTXDim = newXdim/2+1;
	kernel1(imgs, frame.data.xdim, frame.data.ydim, noOfImgs, newXdim, newYdim, filter.data, tmpResult);
// 	******************
//	FIXME normalization has to be done using original img size, i.e frame.data.xdim*frame.data.ydim
//	******************

	MultidimArray<std::complex<double> > V(1, 1, newYdim, newFFTXDim);
	for (size_t i = 0; i < (newFFTXDim*newYdim); i++) {
		V.data[i].real() = tmpResult[i].real() / (frame.data.xdim*frame.data.ydim);
		V.data[i].imag() = tmpResult[i].imag() / (frame.data.xdim*frame.data.ydim);
	}
	Image<double> aaa(newFFTXDim, newYdim, 1, noOfImgs);
	for (size_t i = 0; i < (newFFTXDim*newYdim*noOfImgs); i++) {
		double d = tmpResult[i].real() / (frame.data.xdim*frame.data.ydim);
		if (d < 3) aaa.data[i] = d;
	}
	aaa.write("fftFromGPU.vol");
	std::cout << "normalization done" << std::endl;
	Image<double> yyy (newXdim, newYdim, 1, 1);
	FourierTransformer transformer;
	std::cout << "about to do IFFT" << std::endl;
	transformer.inverseFourierTransform(V, yyy.data);
	std::cout << "IFFT done" << std::endl;
	yyy.write("filteredCroppedInputGPU0.vol");


	// 16785408 X:2049 Y:4096
//	Image<float> tmp(newFFTXDim, newYdim, 1, noOfImgs);
//	for (size_t i = 0; i < (newFFTXDim*newYdim*2); i++) {
////	for (size_t i = 0; i < 8388608L; i++) {
//		float val = result[i].real() / (newYdim*newYdim);
//		if (val < 3) tmp.data[i] = val;
//		else std::cout << "skipping " << val << " at position " << i << std::endl;
//
//	}
//	tmp.write("fftFromGPU" + SSTR(counter) + ".vol");

}

void ProgMovieAlignmentCorrelationGPU::computeShifts(size_t N,
		const Matrix1D<double>& bX, const Matrix1D<double>& bY,
		const Matrix2D<double>& A) {

	float* result;
	kernel3(maxShift, N, tmpResult, newXdim/2+1, newYdim, result);
	std::cout << "kernel3 done" << std::endl;
	size_t framexdim = 4096;
	size_t frameydim = 4096; // FIXME


	size_t newFFTXDim = newXdim/2+1;
	for (int img = 0; img < (N * (N-1)/2); img++) {
//		MultidimArray<std::complex<double> > V(1, 1, newYdim, newFFTXDim);
//		for (size_t i = 0; i < (newFFTXDim*newYdim); i++) {
//			V.data[i].real() = result[i + img*newYdim*newFFTXDim].real() / (framexdim*frameydim);
//			V.data[i].imag() = result[i + img*newYdim*newFFTXDim].imag() / (framexdim*frameydim);
//		}
//		std::cout << "V done" << std::endl;
//		Image<double> aaa(newFFTXDim, newYdim, 1, 1);
//		for (size_t i = 0; i < (newFFTXDim*newYdim); i++) {
//			double d = result[i + img*newYdim*newFFTXDim].real() / (framexdim*frameydim);
//			if (d < 3) aaa.data[i] = d;
//		}
//		aaa.write("correlationGPU" + SSTR(img) + ".vol");
//		std::cout << "correlation done" << std::endl;
//		Image<double> yyy (newXdim, newYdim, 1, 1);
//		FourierTransformer transformer;
//		std::cout << "about to do IFFT" << std::endl;
//		transformer.inverseFourierTransform(V, yyy.data);
//		std::cout << "IFFT done" << std::endl;
//		CenterFFT(yyy.data, true);
//		yyy.write("correlationIFFTGPU" + SSTR(img) + ".vol");
		Image<float>tmp(newXdim, newYdim, 1, 1);
		tmp.data.data = result;
		CenterFFT(tmp.data, true);
		tmp.write("correlationIFFTGPU" + SSTR(img) + ".vol");
	}


	return;
	// FIXME refactor

	int idx = 0;
	MultidimArray<double> Mcorr;
	Mcorr.resizeNoCopy(newYdim, newXdim);
	Mcorr.setXmippOrigin();
	CorrelationAux aux;
	for (size_t i = 0; i < N - 1; ++i) {
		for (size_t j = i + 1; j < N; ++j) {
			bestShift(*frameFourier[i], *frameFourier[j], Mcorr, bX(idx),
					bY(idx), aux, NULL, maxShift);
			if (verbose)
				std::cerr << "Frame " << i + nfirst << " to Frame "
						<< j + nfirst << " -> (" << bX(idx) << "," << bY(idx)
						<< ")\n";
			for (int ij = i; ij < j; ij++)
				A(idx, ij) = 1;

			idx++;
		}
		delete frameFourier[i];
	}
}
