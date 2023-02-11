/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package volvis;

import com.jogamp.opengl.GL;
import com.jogamp.opengl.GL2;
import com.jogamp.opengl.util.texture.Texture;
import com.jogamp.opengl.util.texture.awt.AWTTextureIO;
import gui.RaycastRendererPanel;
import gui.TransferFunction2DEditor;
import gui.TransferFunctionEditor;
import java.awt.image.BufferedImage;
import util.TFChangeListener;
import util.VectorMath;
import volume.GradientVolume;
import volume.Volume;
//add one more package:volume.VoxelGradient
import volume.VoxelGradient;

/**
 *
 * @author michel
 */
public class RaycastRenderer extends Renderer implements TFChangeListener {

    private Volume volume = null;
    private GradientVolume gradients = null;
    RaycastRendererPanel panel;
    TransferFunction tFunc;
    TransferFunctionEditor tfEditor;
    TransferFunction2DEditor tfEditor2D;
    //add some variables
    private int mode=1;
    private double alpha=10;
    private double k_ambient=0.1;
    private double k_diffuse=0.7;
    private double k_specular=0.2;
    private TFColor ambient=new TFColor();
    private TFColor diffuse=new TFColor();
    private boolean is_shading=false;  // checkbox for shading, ture for ticking,false for not ticking
	//variables added.
	public RaycastRenderer() {
        panel = new RaycastRendererPanel(this);
        panel.setSpeedLabel("0");
    }
	//add some functions which will be used in 'RaycastRendererPanel.java'
 	public void setMode(int mode){
        this.mode=mode;
        this.changed();
    }
    public void switchShading(){
     is_shading=(!is_shading);
     this.changed();
    }
    //functions about 'RaycastRendererPanel.java' are added.
    
    
	//add functions of rendering algorithms
    //pack the process of Phong shading into getShadedColor(******).
    private TFColor getShadedColor(double[] coord, TFColor voxelColor, double[] viewVec) {
        if (coord[0]<0 || volume.getDimX()<=coord[0]  
         || coord[1]<0 || volume.getDimY()<=coord[1]
         || coord[2]<0 || volume.getDimZ()<=coord[2]) {
            voxelColor.r=0;
            voxelColor.g=0;
            voxelColor.b=0;
            return voxelColor;
        }
        VoxelGradient VoxGrad=gradients.getGradient(
                (int) Math.floor(coord[0]),
                (int) Math.floor(coord[1]),
                (int) Math.floor(coord[2]));
        double[] N=VoxGrad.getNormal();
        double L_dotproduct_N=Math.max(0, VectorMath.dotproduct(viewVec, N));
        double[] R=VectorMath.subtract((VectorMath.scale(N, 2*VectorMath.dotproduct(N, viewVec))), viewVec);
        double V_dotproduct_R=VectorMath.dotproduct(viewVec, R);
      	V_dotproduct_R=(V_dotproduct_R<=0 || L_dotproduct_N<=0) ? 0 : Math.pow(V_dotproduct_R, alpha); 
        ambient = new TFColor(1,1,1,1);
        diffuse=voxelColor;
        //here comes the formula
        return new TFColor(k_ambient*ambient.r+k_diffuse*diffuse.r*L_dotproduct_N+k_specular*diffuse.r*V_dotproduct_R,
						   k_ambient*ambient.g+k_diffuse*diffuse.g*L_dotproduct_N+k_specular*diffuse.g*V_dotproduct_R,
					 	   k_ambient*ambient.b+k_diffuse*diffuse.b*L_dotproduct_N+k_specular*diffuse.b*V_dotproduct_R,
                		   diffuse.a);
    }
    //shading function added.

    public void setVolume(Volume vol) {
        System.out.println("Assigning volume");
        volume = vol;

        System.out.println("Computing gradients");
        gradients = new GradientVolume(vol);

        // set up image for storing the resulting rendering
        // the image width and height are equal to the length of the volume diagonal
        int imageSize = (int) Math.floor(Math.sqrt(vol.getDimX() * vol.getDimX() + vol.getDimY() * vol.getDimY()
                + vol.getDimZ() * vol.getDimZ()));
        if (imageSize % 2 != 0) {
            imageSize = imageSize + 1;
        }
        image = new BufferedImage(imageSize, imageSize, BufferedImage.TYPE_INT_ARGB);
        // create a standard TF where lowest intensity maps to black, the highest to white, and opacity increases
        // linearly from 0.0 to 1.0 over the intensity range
        tFunc = new TransferFunction(volume.getMinimum(), volume.getMaximum());
        
        // uncomment this to initialize the TF with good starting values for the orange dataset 
        tFunc.setTestFunc();
        
        
        tFunc.addTFChangeListener(this);
        tfEditor = new TransferFunctionEditor(tFunc, volume.getHistogram());
        
        tfEditor2D = new TransferFunction2DEditor(volume, gradients);
        tfEditor2D.addTFChangeListener(this);

        System.out.println("Finished initialization of RaycastRenderer");
    }

    public RaycastRendererPanel getPanel() {
        return panel;
    }

    public TransferFunction2DEditor getTF2DPanel() {
        return tfEditor2D;
    }
    
    public TransferFunctionEditor getTFPanel() {
        return tfEditor;
    }
    //modify getVoxel(***) function, add trilinear(tri-linear) feature to this function
    double getVoxel(double[] coord) {
        	double x = coord[0];
            double y = coord[1];
            double z = coord[2];

            int x_floor = (int) Math.floor(x);
            int x_ceil = (int) Math.min(Math.ceil(x), volume.getDimX() - 1);
            int y_floor = (int) Math.floor(y);
            int y_ceil = (int) Math.min(Math.ceil(y), volume.getDimY() - 1);
            int z_floor = (int) Math.floor(z);
            int z_ceil = (int) Math.min(Math.ceil(z), volume.getDimZ() - 1);

            if (coord[0] < 0 || volume.getDimX() <= coord[0]  
             		|| coord[1] < 0 || volume.getDimY() <= coord[1]
                    || coord[2] < 0 || volume.getDimZ() <= coord[2] 
                    || x_floor < 0 || volume.getDimX() <= x_floor  
                    || y_floor < 0 || volume.getDimY() <= y_floor
                    || z_floor < 0 || volume.getDimZ() <= z_floor
                    || x_ceil < 0 || volume.getDimX() <= x_ceil 
                    || y_ceil < 0 || volume.getDimY() <= y_ceil 
                    || z_ceil < 0 || volume.getDimZ() <= z_ceil  ){
                return 0;
            }
            //start to add trilinear feature
			//get voxels measured and needed            
            double vfff = volume.getVoxel(x_floor, y_floor, z_floor);
            double vffc = volume.getVoxel(x_floor, y_floor, z_ceil);
            double vfcf = volume.getVoxel(x_floor, y_ceil, z_floor);
            double vfcc = volume.getVoxel(x_floor, y_ceil, z_ceil);
            double vcff = volume.getVoxel(x_ceil, y_floor, z_floor);
            double vcfc = volume.getVoxel(x_ceil, y_floor, z_ceil);
            double vccf = volume.getVoxel(x_ceil, y_ceil, z_floor);
            double vccc = volume.getVoxel(x_ceil, y_ceil, z_ceil);

            //estimations on x dimension
            double vff = ((1-((x-x_floor)/(x_ceil-x_floor)))*vfff )+ (((x-x_floor)/(x_ceil-x_floor))*vcff);
            double vfc = ((1-((x-x_floor)/(x_ceil-x_floor)))*vffc )+ (((x-x_floor)/(x_ceil-x_floor))*vcfc);
            double vcf = ((1-((x-x_floor)/(x_ceil-x_floor)))*vfcf )+ (((x-x_floor)/(x_ceil-x_floor))*vccf);
            double vcc = ((1-((x-x_floor)/(x_ceil-x_floor)))*vfcc )+ (((x-x_floor)/(x_ceil-x_floor))*vccc);

            // estimations on y dimension
            double vf = ((1-((y-y_floor)/(y_ceil-y_floor)))*vff )+ (((y-y_floor)/(y_ceil-y_floor))*vcf);
            double vc = ((1-((y-y_floor)/(y_ceil-y_floor)))*vfc )+ (((y-y_floor)/(y_ceil-y_floor))*vcc);

            // get final estimation on z dimension
            double res = ((1-((z-z_floor)/(z_ceil-z_floor)))*vf )+ (((z-z_floor)/(z_ceil-z_floor))*vc);
            //trilinear features added
            return res;
        
    }  
	//function modified.
    void slicer(double[] viewMatrix) {
    	//add a line
        int step_length =interactiveMode?5:1;
        // clear image
        for (int j = 0; j < image.getHeight(); j++) {
            for (int i = 0; i < image.getWidth(); i++) {
                image.setRGB(i, j, 0);
            }
        }

        // vector uVec and vVec define a plane through the origin, 
        // perpendicular to the view vector viewVec
        double[] viewVec = new double[3];
        double[] uVec = new double[3];
        double[] vVec = new double[3];
        VectorMath.setVector(viewVec, viewMatrix[2], viewMatrix[6], viewMatrix[10]);
        VectorMath.setVector(uVec, viewMatrix[0], viewMatrix[4], viewMatrix[8]);
        VectorMath.setVector(vVec, viewMatrix[1], viewMatrix[5], viewMatrix[9]);

        // image is square
        int imageCenter = image.getWidth() / 2;

        double[] pixelCoord = new double[3];
        double[] volumeCenter = new double[3];
        VectorMath.setVector(volumeCenter, volume.getDimX() / 2, volume.getDimY() / 2, volume.getDimZ() / 2);

        // sample on a plane through the origin of the volume data
        double max = volume.getMaximum();
        TFColor voxelColor = new TFColor();

        for (int j = 0; j < image.getHeight(); j=j+step_length) {
        	for (int i = 0; i < image.getWidth(); i=i+step_length) {
                pixelCoord[0]=uVec[0]*(i-imageCenter)+vVec[0]*(j-imageCenter)+volumeCenter[0];
                pixelCoord[1]=uVec[1]*(i-imageCenter)+vVec[1]*(j-imageCenter)+volumeCenter[1];
				pixelCoord[2]=uVec[2]*(i-imageCenter)+vVec[2]*(j-imageCenter)+volumeCenter[2];

                double val = getVoxel(pixelCoord);
                
                // Map the intensity to a grey value by linear scaling
                voxelColor.r = val/max;
                voxelColor.g = voxelColor.r;
                voxelColor.b = voxelColor.r;
                voxelColor.a = val > 0 ? 1.0 : 0.0;  // this makes intensity 0 completely transparent and the rest opaque
                // Alternatively, apply the transfer function to obtain a color
                //voxelColor = tFunc.getColor(val);
                
                
                // BufferedImage expects a pixel color packed as ARGB in an int
                int c_alpha = voxelColor.a <= 1.0 ? (int) Math.floor(voxelColor.a * 255) : 255;
                int c_red = voxelColor.r <= 1.0 ? (int) Math.floor(voxelColor.r * 255) : 255;
                int c_green = voxelColor.g <= 1.0 ? (int) Math.floor(voxelColor.g * 255) : 255;
                int c_blue = voxelColor.b <= 1.0 ? (int) Math.floor(voxelColor.b * 255) : 255;
                int pixelColor = (c_alpha << 24) | (c_red << 16) | (c_green << 8) | c_blue;
                //modify setRGB
                for (int u = 0; u < step_length; u++) {
                    for (int v= 0; v < step_length; v++) {
                        if (v + i < image.getHeight() 
                            && v + i >= 0
                            && u + j < image.getWidth() 
                            && u + j >= 0) {
                            image.setRGB(v + i, u + j,pixelColor);
                        }
                    }
                }
            }
        }

    }
    
      

	void MIP(double[] viewMatrix) {
        int step_length = interactiveMode?5:1;
        // clear image
        for (int j = 0; j < image.getHeight(); j++) {
            for (int i = 0; i < image.getWidth(); i++) {
                image.setRGB(i, j, 0);
            }
        }
        double[] uVec = new double[3];
        double[] vVec = new double[3];
        double[] viewVec = new double[3];     
        // vector uVec and vVec define a plane through the origin, 
        // perpendicular to the view vector viewVec
        VectorMath.setVector(viewVec, viewMatrix[2], viewMatrix[6], viewMatrix[10]);
        VectorMath.setVector(uVec, viewMatrix[0], viewMatrix[4], viewMatrix[8]);
        VectorMath.setVector(vVec, viewMatrix[1], viewMatrix[5], viewMatrix[9]);
		//image is square
        int imageCenter = image.getWidth() / 2;
        double[] pixelCoord = new double[3];
        double[] volumeCenter = new double[3];
        VectorMath.setVector(volumeCenter, volume.getDimX() / 2, volume.getDimY() / 2, volume.getDimZ() / 2);
		//sample on a plane through the origin of the volume data
        double maxOfAll = volume.getMaximum();
        TFColor voxelColor = new TFColor();
        // traverse along the ray casted from the pixel, whose cooridination in
        // u,v system is (i - imageCenter,j - imageCenter). During the traversing process, 
        // the coordination in u,v,view 3d system is 
        // (i - imageCenter,j - imageCenter,step - imageCenter)
        // the origin point for the u,v,view 3d system is at the center of the voxel.
        // the origin point for the voxel 3d system is at a cube vertex of the voxel.
        // the x,y,z vectors for the voxel 3d system are cube sides of the voxel.
		for (int j = 0; j < image.getHeight(); j+= step_length) {
        	for (int i = 0; i < image.getWidth(); i+= step_length) {
                int diagonal=(int)Math.sqrt(volume.getDimX()*volume.getDimX()+volume.getDimY()*volume.getDimY()+volume.getDimZ()*volume.getDimZ()); 
                double maxOnRay = 0;
                for (int step = 0; step < diagonal - 1; step+= step_length) {
					pixelCoord[0]=uVec[0]*(i-imageCenter)+vVec[0]*(j-imageCenter)+(step-imageCenter)*viewVec[0]+volumeCenter[0];
					pixelCoord[1]=uVec[1]*(i-imageCenter)+vVec[1]*(j-imageCenter)+(step-imageCenter)*viewVec[1]+volumeCenter[1];
                    pixelCoord[2]=uVec[2]*(i-imageCenter)+vVec[2]*(j-imageCenter)+(step-imageCenter)*viewVec[2]+volumeCenter[2];
					//Do not need to consider the validity of the pixelCoord, because it will
					// be considered in this.getVoxel(***) function. If invalid, return 0.
					double val = getVoxel(pixelCoord);
                    if (val > maxOnRay) {
                        maxOnRay = val;
                    }
                    if (val / maxOfAll > 0.93) {
                        break;
                    }
                }
               
                
                // Map the intensity to a grey value by linear scaling
                voxelColor.r = maxOnRay/maxOfAll;
                voxelColor.g = voxelColor.r;
                voxelColor.b = voxelColor.r;
                voxelColor.a = maxOnRay > 0 ? 1.0 : 0.0;  // this makes intensity 0 completely transparent and the rest opaque
                // Alternatively, apply the transfer function to obtain a color
                //voxelColor = tFunc.getColor(val);
                
                
                // BufferedImage expects a pixel color packed as ARGB in an int
                int c_alpha = voxelColor.a <= 1.0 ? (int) Math.floor(voxelColor.a * 255) : 255;
                int c_red = voxelColor.r <= 1.0 ? (int) Math.floor(voxelColor.r * 255) : 255;
                int c_green = voxelColor.g <= 1.0 ? (int) Math.floor(voxelColor.g * 255) : 255;
                int c_blue = voxelColor.b <= 1.0 ? (int) Math.floor(voxelColor.b * 255) : 255;
                int pixelColor = (c_alpha << 24) | (c_red << 16) | (c_green << 8) | c_blue;
                for (int u = 0; u < step_length; u++) {
                    for (int v= 0; v < step_length; v++) {
                        if (v + i < image.getHeight() 
                            && v + i >= 0
                            && u + j < image.getWidth() 
                            && u + j >= 0) {
                            image.setRGB(v + i, u + j,pixelColor);
                        }
                    }
                }
            }
        }

    }
    private void Composite(double[] viewMatrix) {
        int step_length =interactiveMode?5:1;
        // vector uVec and vVec define a plane through the origin, 
        // perpendicular to the view vector viewVec
        double[] uVec = new double[3];
        double[] vVec = new double[3];
        double[] viewVec = new double[3];
        VectorMath.setVector(uVec, viewMatrix[0], viewMatrix[4], viewMatrix[8]);
        VectorMath.setVector(vVec, viewMatrix[1], viewMatrix[5], viewMatrix[9]);
        VectorMath.setVector(viewVec, viewMatrix[2], viewMatrix[6], viewMatrix[10]);
        VectorMath.normalize(viewVec);
        // image is square
        int imageCenter = image.getWidth() / 2;

        
        double[] volumeCenter = new double[3];
        VectorMath.setVector(volumeCenter, volume.getDimX() / 2, volume.getDimY() / 2, volume.getDimZ() / 2);

        // sample on a plane through the origin of the volume data
        TFColor voxelColorHere, ColorOnTheRay = new TFColor();
    	double[] pixelCoord = new double[3];
        for (int j = 0; j < image.getHeight(); j=j+step_length) {
            for (int i = 0; i < image.getWidth(); i=i+step_length) {
                ColorOnTheRay.r=0;
                ColorOnTheRay.g=0;
                ColorOnTheRay.b=0;
                ColorOnTheRay.a=1;
                int diagonal = (int) Math.sqrt(volume.getDimX() * volume.getDimX() + volume.getDimY() * volume.getDimY() + volume.getDimZ() * volume.getDimZ()); 
                for (int step =0;step < diagonal - 1;step+=step_length) {
	                pixelCoord[0]=uVec[0]*(i-imageCenter)+vVec[0] * (j - imageCenter) + viewVec[0]*(step - imageCenter ) + volumeCenter[0];
	                pixelCoord[1] = uVec[1] * (i - imageCenter)
	                        + vVec[1] * (j - imageCenter) 
	                         + viewVec[1]*(step - imageCenter) + volumeCenter[1];
	                pixelCoord[2] = uVec[2] * (i - imageCenter) 
	                        + vVec[2] * (j - imageCenter)
	                         + viewVec[2]*(step - imageCenter)+volumeCenter[2];
					double val = getVoxel(pixelCoord);

	                // apply the transfer function to obtain a color
	                voxelColorHere = tFunc.getColor(val);
					//phong shading
                    if (is_shading) {
                        voxelColorHere = getShadedColor(pixelCoord, voxelColorHere, viewVec);
                    }
                    //compositing
                    ColorOnTheRay.r = voxelColorHere.a * voxelColorHere.r + (1 - voxelColorHere.a) * ColorOnTheRay.r;
                    ColorOnTheRay.g = voxelColorHere.a * voxelColorHere.g + (1 - voxelColorHere.a) * ColorOnTheRay.g;
                    ColorOnTheRay.b = voxelColorHere.a * voxelColorHere.b + (1 - voxelColorHere.a) * ColorOnTheRay.b;
                }
 				//BufferedImage expects a pixel color packed as ARGB in an int
                int c_alpha = ColorOnTheRay.a <= 1.0 ? (int) Math.floor(ColorOnTheRay.a * 255) : 255;
                int c_red = ColorOnTheRay.r <= 1.0 ? (int) Math.floor(ColorOnTheRay.r * 255) : 255;
                int c_green = ColorOnTheRay.g <= 1.0 ? (int) Math.floor(ColorOnTheRay.g * 255) : 255;
                int c_blue = ColorOnTheRay.b <= 1.0 ? (int) Math.floor(ColorOnTheRay.b * 255) : 255;
                int pixelColor = (c_alpha << 24) | (c_red << 16) | (c_green << 8) | c_blue;
                for (int u = 0; u <step_length;u++) {
                    for (int v= 0; v <step_length; v++) {
                        if (v + i < image.getHeight() 
                            && v + i >= 0
                            && u + j < image.getWidth() 
                            && u + j >= 0) {
                            image.setRGB(v + i, u + j,pixelColor);
                        }
                    }
                }
            }
        }
    }
     void Transfer2DFunction (double[] viewMatrix){
        int step_length =interactiveMode?5:1;
        // clear image
        for (int j = 0; j < image.getHeight(); j++) {
            for (int i = 0; i < image.getWidth(); i++) {
                image.setRGB(i, j, 0);
            }
        }
		// vector uVec and vVec define a plane through the origin, 
        // perpendicular to the view vector viewVec
        double[] viewVec = new double[3];
        double[] uVec = new double[3];
        double[] vVec = new double[3];  
        VectorMath.setVector(viewVec, viewMatrix[2], viewMatrix[6], viewMatrix[10]);
        VectorMath.setVector(uVec, viewMatrix[0], viewMatrix[4], viewMatrix[8]);
        VectorMath.setVector(vVec, viewMatrix[1], viewMatrix[5], viewMatrix[9]);
		// image is square
        int imageCenter = image.getWidth() / 2;
        double[] pixelCoord = new double[3];
        double[] volumeCenter = new double[3];
        VectorMath.setVector(volumeCenter, volume.getDimX() / 2, volume.getDimY() / 2, volume.getDimZ() / 2);
		// sample on a plane through the origin of the volume data
        double max = volume.getMaximum();
        // Get the maximum length that is diagonal
        int diagonal = (int) Math.sqrt(volume.getDimX() * volume.getDimX() + volume.getDimY() * volume.getDimY() + volume.getDimZ() * volume.getDimZ());
        // For each pixel of the image
        for (int j = 0; j < image.getHeight(); j+=step_length) {
            for (int i = 0; i < image.getWidth(); i+=step_length) {
                // First, set a color variable in which we can put all the colors together
                TFColor ColorOnTheRay = new TFColor(0,0,0,1);
				for(int step = 0; step < diagonal - 1; step+=step_length) {
                    pixelCoord[0] = uVec[0] * (i - imageCenter) + vVec[0] * (j - imageCenter) + viewVec[0] * (step - imageCenter) + volumeCenter[0];
                    pixelCoord[1] = uVec[1] * (i - imageCenter) + vVec[1] * (j - imageCenter) + viewVec[1] * (step - imageCenter) + volumeCenter[1];
                    pixelCoord[2] = uVec[2] * (i - imageCenter) + vVec[2] * (j - imageCenter) + viewVec[2] * (step - imageCenter) + volumeCenter[2];
                  	int val = (int) getVoxel(pixelCoord);
                    if (val == 0) {
                        continue;
                    }
                    TransferFunction2DEditor.TriangleWidget triagWid = this.tfEditor2D.triangleWidget;
                    float baseIntensity = triagWid.baseIntensity;
                    TFColor voxelColorHere = triagWid.color;
                    double radius = triagWid.radius;//constant,0.2, at line 50 and line 328 of TransferFunction2DEditor.java
                    double opacity = 0;
                    VoxelGradient gradient = gradients.getGradient((int) pixelCoord[0], (int) pixelCoord[1], (int) pixelCoord[2]);
                    if (val == baseIntensity && gradient.mag == 0){
                        opacity = 1;
                    } else if (gradient.mag > 0 && ((val - (radius * gradient.mag)) <= baseIntensity && baseIntensity <= (val + (radius * gradient.mag)))) {
                            opacity = voxelColorHere.a * (1 - (1 / radius) * Math.abs(((baseIntensity - val) / gradient.mag)));
                        } else {
                        opacity = 0;
                    }
                    if (is_shading) {
                        voxelColorHere = getShadedColor(pixelCoord, voxelColorHere,viewVec);
                    } 
                    ColorOnTheRay.r = (voxelColorHere.r * opacity ) + (ColorOnTheRay.r * (1 - opacity));
                    ColorOnTheRay.g = (voxelColorHere.g * opacity ) + (ColorOnTheRay.g * (1 - opacity));
                    ColorOnTheRay.b = (voxelColorHere.b * opacity ) + (ColorOnTheRay.b * (1 - opacity)); 
                  }
                int c_alpha = ColorOnTheRay.a <= 1.0 ? (int) Math.floor(ColorOnTheRay.a * 255) : 255;
                int c_red = ColorOnTheRay.r <= 1.0 ? (int) Math.floor(ColorOnTheRay.r * 255) : 255;
                int c_green = ColorOnTheRay.g <= 1.0 ? (int) Math.floor(ColorOnTheRay.g * 255) : 255;
                int c_blue = ColorOnTheRay.b <= 1.0 ? (int) Math.floor(ColorOnTheRay.b * 255) : 255;
                int pixelColor = (c_alpha << 24) | (c_red << 16) | (c_green << 8) | c_blue;
                image.setRGB(i, j, pixelColor);
                for (int u = 0; u <step_length; u++) {
                    for (int v= 0; v <step_length; v++) {
                        if (v + i < image.getHeight() 
                            && v + i >= 0
                            && u + j < image.getWidth() 
                            && u + j >= 0) {
                            image.setRGB(v + i, u + j,pixelColor);
                        }
                    }
                }
            }
        }
     }
  
    private void drawBoundingBox(GL2 gl) {
        gl.glPushAttrib(GL2.GL_CURRENT_BIT);
        gl.glDisable(GL2.GL_LIGHTING);
        gl.glColor4d(1.0, 1.0, 1.0, 1.0);
        gl.glLineWidth(1.5f);
        gl.glEnable(GL.GL_LINE_SMOOTH);
        gl.glHint(GL.GL_LINE_SMOOTH_HINT, GL.GL_NICEST);
        gl.glEnable(GL.GL_BLEND);
        gl.glBlendFunc(GL.GL_SRC_ALPHA, GL.GL_ONE_MINUS_SRC_ALPHA);

        gl.glBegin(GL.GL_LINE_LOOP);
        gl.glVertex3d(-volume.getDimX() / 2.0, -volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(-volume.getDimX() / 2.0, volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, -volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glEnd();

        gl.glBegin(GL.GL_LINE_LOOP);
        gl.glVertex3d(-volume.getDimX() / 2.0, -volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glVertex3d(-volume.getDimX() / 2.0, volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, -volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glEnd();

        gl.glBegin(GL.GL_LINE_LOOP);
        gl.glVertex3d(volume.getDimX() / 2.0, -volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, -volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glEnd();

        gl.glBegin(GL.GL_LINE_LOOP);
        gl.glVertex3d(-volume.getDimX() / 2.0, -volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glVertex3d(-volume.getDimX() / 2.0, -volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(-volume.getDimX() / 2.0, volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(-volume.getDimX() / 2.0, volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glEnd();

        gl.glBegin(GL.GL_LINE_LOOP);
        gl.glVertex3d(-volume.getDimX() / 2.0, volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glVertex3d(-volume.getDimX() / 2.0, volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glEnd();

        gl.glBegin(GL.GL_LINE_LOOP);
        gl.glVertex3d(-volume.getDimX() / 2.0, -volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glVertex3d(-volume.getDimX() / 2.0, -volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, -volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, -volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glEnd();

        gl.glDisable(GL.GL_LINE_SMOOTH);
        gl.glDisable(GL.GL_BLEND);
        gl.glEnable(GL2.GL_LIGHTING);
        gl.glPopAttrib();

    }

    @Override
    public void visualize(GL2 gl) {


        if (volume == null) {
            return;
        }

        drawBoundingBox(gl);

        gl.glGetDoublev(GL2.GL_MODELVIEW_MATRIX, viewMatrix, 0);//??

        long startTime = System.currentTimeMillis();
            
            switch (mode){
            case 1:
                slicer(viewMatrix);
                break;
            case 2:
                MIP(viewMatrix);
                break;
            case 3:
                Composite(viewMatrix);
                break;
            case 4:
                Transfer2DFunction(viewMatrix);
                break;
            default:
                slicer(viewMatrix);
        }    
        long endTime = System.currentTimeMillis();
        double runningTime = endTime - startTime;
        panel.setSpeedLabel(Double.toString(runningTime));

        Texture texture = AWTTextureIO.newTexture(gl.getGLProfile(), image, false);

        gl.glPushAttrib(GL2.GL_LIGHTING_BIT);
        gl.glDisable(GL2.GL_LIGHTING);
        gl.glEnable(GL.GL_BLEND);
        gl.glBlendFunc(GL.GL_SRC_ALPHA, GL.GL_ONE_MINUS_SRC_ALPHA);

        // draw rendered image as a billboard texture
        texture.enable(gl);
        texture.bind(gl);
        double halfWidth = image.getWidth() / 2.0;
        gl.glPushMatrix();
        gl.glLoadIdentity();
        gl.glBegin(GL2.GL_QUADS);
        gl.glColor4f(1.0f, 1.0f, 1.0f, 1.0f);
        gl.glTexCoord2d(0.0, 0.0);
        gl.glVertex3d(-halfWidth, -halfWidth, 0.0);
        gl.glTexCoord2d(0.0, 1.0);
        gl.glVertex3d(-halfWidth, halfWidth, 0.0);
        gl.glTexCoord2d(1.0, 1.0);
        gl.glVertex3d(halfWidth, halfWidth, 0.0);
        gl.glTexCoord2d(1.0, 0.0);
        gl.glVertex3d(halfWidth, -halfWidth, 0.0);
        gl.glEnd();
        texture.disable(gl);
        texture.destroy(gl);
        gl.glPopMatrix();

        gl.glPopAttrib();


        if (gl.glGetError() > 0) {
            System.out.println("some OpenGL error: " + gl.glGetError());
        }

    }
    private BufferedImage image;
    private double[] viewMatrix = new double[4 * 4];
    
    @Override
    public void changed() {
        for (int i=0; i < listeners.size(); i++) {
            listeners.get(i).changed();
        }
    }
}
