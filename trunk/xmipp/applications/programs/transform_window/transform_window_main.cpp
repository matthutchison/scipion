/***************************************************************************
 *
 * Authors:     Carlos Oscar S. Sorzano (coss@cnb.csic.es)
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

#include <data/progs.h>
#include <data/args.h>
#include <data/image_generic.h>

class ProgWindow: public XmippMetadataProgram
{
public:
    typedef enum {CORNERMODE, SIZEMODE, CROPMODE} WindowMode;

    bool physical_coords;
    int sizeX, sizeY, sizeZ;
    int cropX, cropY, cropZ;
    int x0, y0, z0;
    int xF, yF, zF;
    double padValue;
    std::string padType;
    WindowMode mode;
    Image<char>                IChar;
    Image<unsigned char>       IUChar;
    Image<short int>           IShort;
    Image<unsigned short int>  IUShort;
    Image<float>               IFloat;

    void defineParams()
    {
        each_image_produces_an_output = true;
        XmippMetadataProgram::defineParams();
        addUsageLine("This program takes a region from the input images");
        addUsageLine(" ");
        addUsageLine("But... what is windowing? Let's suppose we have an image 32x32,");
        addUsageLine("and we window it to 16x16, then the central part which contains");
        addUsageLine("this 16x16 square is returned, but if we window to 64x64 then the");
        addUsageLine("image is kept in the center of the final image and it is padded");
        addUsageLine("with 0's until the new size is reached (the padding value may be");
        addUsageLine("modified using --pad_value, --corner_pad_value or --average_pad_value).");
        addUsageLine(" ");
        addUsageLine("With this idea in mind of padding or cutting, if you define two");
        addUsageLine("logical corners in the input image (the most negative -the top");
        addUsageLine("left in an image- and the most positive -bottom right in an image-)");
        addUsageLine("you can store the new image or volume generated by windowing with");
        addUsageLine("a square box.");
        addParamsLine("  --corners <...>                        : Windows corners (by default indexes are logical)");
        addParamsLine("                                         : 2D: <x0> <y0> <xF> <yF>");
        addParamsLine("                                         : 3D: <x0> <y0> <z0> <xF> <yF> <zF>");
        addParamsLine("  or --size <sizeX> <sizeY=0> <sizeZ=0>  : Output dimension. The volume is windowed");
        addParamsLine("                                         : (expanded or cutted) in all directions such");
        addParamsLine("                                         : that the origin of the volume remains the");
        addParamsLine("                                         : same. If the Y and Z dimensions are not");
        addParamsLine("                                         : specified they are assumed to be the same as");
        addParamsLine("                                         : the X dimension.");
        addParamsLine("  or --crop <sizeX> <sizeY=0> <sizeZ=0>  : Crop this amount of pixels in each direction.");
        addParamsLine("                                         : Half of the pixels will be cropped from the left");
        addParamsLine("                                         : and the other half from the right");
        addParamsLine("                                         : if only one is given, the other two");
        addParamsLine("                                         : are supposed to be the same");
        addParamsLine("  [--physical]                           : use physical instead of logical coordinates");
        addParamsLine("    requires --corners;");
        addParamsLine("  [--pad <padtype=value>]                : value used for padding");
        addParamsLine("   where <padtype>");
        addParamsLine("         value <v=0>                     : use this value for padding");
        addParamsLine("         corner                          : use the top-left corner for padding");
        addParamsLine("         avg                             : use the image average for padding");
        addExampleLine("Window a single image to 16x16, overwriting input image",false);
        addExampleLine("xmipp_transform_window -i g0ta0001.xmp --size 16");
        addExampleLine("Window a single volume to 32x64x64",false);
        addExampleLine("xmipp_transform_window -i g0ta.vol --size 64 64 32");
        addExampleLine("The same using logical indexes",false);
        addExampleLine("xmipp_transform_window -i g0ta.vol --corners -32 -32 -16 31 31 15");
        addExampleLine("Note that r0 and rF are not symmetric because the volume is of an",false);
        addExampleLine("even size, if we wanted to get a 33x65x65 volume, the right indexes",false);
        addExampleLine("would be",false);
        addExampleLine("xmipp_transform_window -i g0ta.vol --corners -32 -32 -16 32 32 16");
        addExampleLine("Reduce the volume by 10 pixels on each direction",false);
        addExampleLine("xmipp_transform_window -i g0ta.vol --crop 10");
        addExampleLine("Enlarge the volume by 10 pixels on each direction (negative crop)",false);
        addExampleLine("xmipp_transform_window -i g0ta.vol --crop -10");
    }

    void readParams()
    {
        XmippMetadataProgram::readParams();
        padType=getParam("--pad");
        if (padType == "value")
            padValue=getDoubleParam("--pad","value");
        if (checkParam("--corners"))
        {
            int nparams=getCountParam("--corners");
            if (nparams==4 || nparams==6)
            {
                x0=getIntParam("--corners",0);
                y0=getIntParam("--corners",1);
                if (nparams==4)
                {
                    xF=getIntParam("--corners",2);
                    yF=getIntParam("--corners",3);
                    zF=z0=0;
                }
                else
                {
                    z0=getIntParam("--corners",2);
                    xF=getIntParam("--corners",3);
                    yF=getIntParam("--corners",4);
                    zF=getIntParam("--corners",5);
                }
            }
            else
                REPORT_ERROR(ERR_ARG_INCORRECT,"Incorrect number of arguments after --corners");
            physical_coords = checkParam("--physical");
            mode=CORNERMODE;
            // TODO Chequear el número de parámetros
        }
        else if (checkParam("--size"))
        {
            sizeX=getIntParam("--size",0);
            sizeY=getIntParam("--size",1);
            sizeZ=getIntParam("--size",2);
            if (sizeY==0)
                sizeY=sizeX;
            if (sizeZ==0)
                sizeZ=sizeX;
            x0 = FIRST_XMIPP_INDEX(sizeX);
            y0 = FIRST_XMIPP_INDEX(sizeY);
            z0 = FIRST_XMIPP_INDEX(sizeZ);
            xF = LAST_XMIPP_INDEX(sizeX);
            yF = LAST_XMIPP_INDEX(sizeY);
            zF = LAST_XMIPP_INDEX(sizeZ);
            mode=SIZEMODE;
            physical_coords = false;
        }
        else if (checkParam("--crop"))
        {
            cropX=getIntParam("--crop",0);
            cropY=getIntParam("--crop",1);
            cropZ=getIntParam("--crop",2);
            if (cropY==0)
                cropY=cropX;
            if (cropZ==0)
                cropZ=cropX;
            mode=CROPMODE;
            physical_coords = false;
        }
    }

    void show()
    {
        if (verbose==0)
            return;
        XmippMetadataProgram::show();
        switch (mode)
        {
        case SIZEMODE:
            std::cout << "New size: (XxYxZ)=" << sizeX << "x" << sizeY << "x"
            << sizeZ << std::endl;
            break;
        case CROPMODE:
            std::cout << "Crop: (XxYxZ)=" << cropX << "x" << cropY << "x"
            << cropZ << std::endl;
            break;
        case CORNERMODE:
            std::cout << "New window: from (z0,y0,x0)=(" << z0 << ","
            << y0 << "," << x0 << ") to (zF,yF,xF)=(" << zF << "," << yF
            << "," << xF << ")\n"
            << "Physical: " << physical_coords << std::endl;
        }
    }

    void processImage(const FileName &fnImg, const FileName &fnImgOut, size_t objId)
    {
        ImageGeneric Iin;
        bool createTempFile=fnImg==fnImgOut;
        if (single_image)
            Iin.read(fnImg,true,-1,true);
        else
            Iin.readApplyGeo(fnImg);

        double init_value(padValue);
        if (padType=="avg")
            init_value=Iin().computeAvg();
        else if (padType=="corner")
            init_value=Iin()(0,0,0,0);

        Iin().setXmippOrigin();
        if (ZSIZE(Iin()())==1)
            zF=z0=0;

        std::string oext=fnImgOut.getExtension();
        DataType dataType=Iin.getDatatype();
        if (oext=="spi" || oext=="xmp" || oext=="vol" || oext=="stk")
            dataType=Float;
        ImageGeneric result(dataType);
        if (mode==CROPMODE)
        {
            int xl=cropX/2;
            int xr=cropX-xl;
            int yl=cropY/2;
            int yr=cropY-yl;
            int zl=cropZ/2;
            int zr=cropZ-zl;
            if (ZSIZE(Iin()())==1)
            {
                zl=zr=0;
            }
            //call to a generic 4D function;
            result.newMappedFile(FINISHINGX(Iin()())-xr-(STARTINGX(Iin()())+xl)+1,
                                 FINISHINGY(Iin()())-yr-(STARTINGY(Iin()())+yl)+1,
                                 FINISHINGZ(Iin()())-zr-(STARTINGZ(Iin()())+zl)+1,
                                 1,
                                 fnImgOut,createTempFile);
            Iin().window(result(),
                         STARTINGZ(Iin()())+zl, STARTINGY(Iin()())+yl,  STARTINGX(Iin()())+xl,
                         FINISHINGZ(Iin()())-zr,FINISHINGY(Iin()())-yr, FINISHINGX(Iin()())-xr);
        }
        else
            if (!physical_coords)
            {
                result.newMappedFile(xF-x0+1,
                                     yF-y0+1,
                                     zF-z0+1,
                                     1,
                                     fnImgOut,
                                     createTempFile);
                Iin().window(result(),z0, y0, x0, zF, yF,xF, init_value);
            }
            else
            {
                result.newMappedFile(xF-x0+1,
                                     yF-y0+1,
                                     zF-z0+1,
                                     1,
                                     fnImgOut,
                                     createTempFile);
                Iin().window(result(),
                             STARTINGZ(Iin()()) + z0,
                             STARTINGY(Iin()()) + y0,
                             STARTINGX(Iin()()) + x0,
                             STARTINGZ(Iin()()) + zF,
                             STARTINGY(Iin()()) + yF,
                             STARTINGX(Iin()()) + xF,
                             init_value);
            }
        Iin.clear(); // Close the input file
        result.write(fnImgOut);
    }
};

int main(int argc, char **argv)
{
    ProgWindow prm;
    prm.read(argc, argv);
    prm.tryRun();
}
