
import ij.IJ;
import ij.ImagePlus;
import ini.trakem2.display.Patch;
import java.awt.Point;
import java.awt.Rectangle;
import java.awt.geom.AffineTransform;
import java.awt.image.BufferedImage;
import java.io.File;
import java.util.List;
import xmipp.Filename;
import xmipp.ImageDouble;

/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
/**
 *
 * @author Juanjo Vega
 */
public class ImagesBuilder {
    /*
    public static ImagePlus buildStack(List<Patch> patches, Rectangle box) {
    ImageStack is = new ImageStack(box.width, box.height);
    
    for (int i = 0; i < patches.size(); i++) {
    Patch patch = patches.get(i);
    
    try {
    float slice[] = buildSlice(patch, box);
    
    is.addSlice(patch.getImagePlus().getTitle(), slice);
    } catch (Exception e) {
    IJ.log("Exception: " + e.getMessage());
    }
    }
    
    ImagePlus ip = new ImagePlus("", is);
    
    return ip;
    }*/

    public static boolean saveResult(List<Patch> patches, Rectangle box, String path) {

        try {
            // Translation to the top left corner.
            AffineTransform translation = new AffineTransform();
            translation.translate(box.x, box.y);

            // Creates a new ImagePlus big enough.
            int w = box.width;
            int h = box.height;

            float result[] = new float[w * h];//buildSlice(patches.get(0), box);

            // Puts images, pixels, translations in arrays, a much convenient way to use them later.
            Patch patch0 = patches.get(0);
            Patch patch1 = patches.get(1);

            ImagePlus images[] = new ImagePlus[]{
                patch0.getImagePlus(),
                patch1.getImagePlus()
            };

            float pixels[][] = new float[][]{
                (float[]) images[0].getProcessor().convertToFloat().getPixels(),
                (float[]) images[1].getProcessor().convertToFloat().getPixels()
            };

            int sourceW[] = new int[]{
                images[0].getWidth(),
                images[1].getWidth()
            };

            int sourceH[] = new int[]{
                images[0].getHeight(),
                images[1].getHeight()
            };

            Rectangle rectangles[] = new Rectangle[]{
                new Rectangle(sourceW[0], sourceH[0]),
                new Rectangle(sourceW[1], sourceH[1])
            };

            // Transforms concatenated with translation to the box top left corner.
            AffineTransform transforms[] = new AffineTransform[]{
                patch0.getAffineTransform().createInverse(),
                patch1.getAffineTransform().createInverse()
            };

            transforms[0].concatenate(translation);
            transforms[1].concatenate(translation);

            for (int j = 0; j < h; j++) {
                for (int i = 0; i < w; i++) {
                    int imageIndex = 0;
                    boolean done = false;
                    Point currentPoint = new Point(i, j);

                    do {
                        Point sourcePoint = new Point();

                        transforms[imageIndex].transform(currentPoint, sourcePoint);

                        int index = getRectangleIndexContaining(rectangles, sourcePoint);

                        //if (rectangles[imageIndex].contains(sourcePoint)) {
                        if (index == imageIndex) {
                            //assign value...
                            result[j * w + i] = pixels[imageIndex][sourcePoint.y * sourceW[imageIndex] + sourcePoint.x];

                            done = true;
                        } else {
                            result[j * w + i] = Float.MIN_VALUE;
                        }

                        imageIndex++;
                    } while (!done && imageIndex < pixels.length);
                }
            }

            // Writes slice to disk.            
            System.err.println(" >>> Saving result: " + path);
            ImageDouble image = new ImageDouble();
            image.setData(box.width, box.height, 1, toDouble(result));
            image.write(path);
        } catch (Exception ex) {
            ex.printStackTrace(System.err);
            IJ.error("Exception: " + ex.getMessage());

            return false;
        }

        return true;
    }

    static int getRectangleIndexContaining(Rectangle rectangles[], Point point) {
        int contains = -1;
        int i = 0;

        for (; i < rectangles.length; i++) {
            Rectangle rectangle = rectangles[i];
            if (rectangle.contains(point)) {
                if (contains >= 0) {    // More than one rectangle contains the point.
                    return -1;
                }

                contains = i;
            }
        }

        return contains;
    }

    public static boolean saveStack(List<Patch> patches, Rectangle box, String path) {
        File f = new File(path);

        // Removes previous (to avoid adding to a stack with different size images)
        if (f.exists()) {
            f.delete();
        }

        for (int i = 1; i <= patches.size(); i++) {
            Patch patch = patches.get(i - 1);

            try {
                float slice[] = buildSlice(patch, box);

                ImageDouble image = new ImageDouble();
                image.setData(box.width, box.height, 1, toDouble(slice));
//System.err.println(" -> X: " + image.getXsize()+ " --> "+box.width);
//System.err.println(" -> Y: " + image.getYsize()+ " --> "+box.height);
//System.err.println(" -> Z: " + image.getNsize()+ " / "+patches.size());

                // Writes slice to disk.
                System.err.println(" >>> Saving slice: " + i + Filename.SEPARATOR + path);
                image.write(i + Filename.SEPARATOR + path);
            } catch (Exception ex) {
                ex.printStackTrace(System.err);
                IJ.error("Exception: " + ex.getMessage());
                return false;
            }
        }

        return true;
    }

    static double[] toDouble(float[] input) {
        double output[] = null;

        if (input != null) {
            output = new double[input.length];
            for (int i = 0; i < input.length; i++) {
                output[i] = input[i];
            }
        }

        return output;
    }

    static float[] buildSlice(Patch patch, Rectangle box) throws Exception {
        ImagePlus ip = patch.getImagePlus();
        AffineTransform translation = new AffineTransform();
        translation.translate(box.x, box.y);

        AffineTransform impT = patch.getAffineTransform();
        impT.invert(); // Inverts transformation.

        AffineTransform T = new AffineTransform(impT);
        T.concatenate(translation);

        float pixels[] = (float[]) ip.getProcessor().convertToFloat().getPixels();
        Rectangle r = new Rectangle(ip.getWidth(), ip.getHeight());

        BufferedImage bimage = new BufferedImage(box.width, box.height,
                BufferedImage.TYPE_INT_ARGB);
        ImagePlus imp = new ImagePlus("tmp", bimage);

        int w = box.width;
        int h = box.height;
        float slice[] = (float[]) imp.getProcessor().convertToFloat().getPixels();

        Point p = new Point();

        for (int j = 0; j < h; j++) {   // Rows
            p.y = j;
            for (int i = 0; i < w; i++) {   // Columns
                p.x = i;

                Point p_ = new Point();
                T.transform(p, p_);

                // If point is inside image...
                if (r.contains(p_)) {

                    // Bilinear Interpolation.
                    float value = Interpolator.bilinear(pixels, ip.getWidth(), ip.getHeight(), p_.getX(), p_.getY());

                    // Store point.
                    slice[j * w + i] = value;
                }
            }
        }

        return slice;
    }
    /*
    public static ImagePlus sumStack(ImagePlus imp) {
    int w = imp.getWidth();
    int h = imp.getHeight();
    float pixels[] = new float[w * h];
    
    for (int i = 0; i < imp.getStackSize(); i++) {
    float current[] = (float[]) imp.getStack().getProcessor(i + 1).getPixels();
    
    pixels = sumArrays(pixels, current);
    }
    
    FloatProcessor ip = new FloatProcessor(w, h);
    ip.setPixels(pixels);
    
    return new ImagePlus("", ip);
    }
    
    static float[] sumArrays(float[] a, float[] b) {
    float result[] = new float[a.length];
    
    for (int i = 0; i < result.length; i++) {
    result[i] = a[i] + b[i];
    }
    
    return result;
    }*/
}
