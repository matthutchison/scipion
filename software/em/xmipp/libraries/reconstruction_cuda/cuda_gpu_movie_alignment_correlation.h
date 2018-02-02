

void kernel1(float* imgs, size_t oldX, size_t oldY, int noOfImages, size_t newX, size_t newY,
		float* filter,
		std::complex<float>*& result);
//		float*& result);

void kernel3(float maxShift, size_t noOfImgs,
		const std::complex<float>* imgs, size_t fftXdim, size_t fftYdim,
		float*& result);

