package model;

import ij.ImagePlus;
import ij.gui.Roi;
import ij.process.ImageProcessor;

import java.awt.Rectangle;
import java.util.List;

import javax.swing.ImageIcon;

public class Particle {
	
	private int x;
	private int y;
	private boolean dragged;
	private Family family;
	private Micrograph micrograph;
	
	
	public Particle(int x, int y, Family family, Micrograph micrograph)
	{
		this.x = x;
		this.y = y;
		this.family = family;
		this.micrograph = micrograph;
	}


	public Micrograph getMicrograph() {
		return micrograph;
	}


	public void setMicrograph(Micrograph micrograph) {
		this.micrograph = micrograph;
	}


	public Family getFamily() {
		return family;
	}

	public void setFamily(Family family) {
		this.family = family;
	}

	public int getX() {
		return x;
	}

	public void setX(int x) {
		this.x = x;
	}

	public int getY() {
		return y;
	}

	public void setY(int y) {
		this.y = y;
	}
	
	public boolean contains(int size, int x2, int y2 )
	{
		int radius = size/2;
			if(x2 < x - radius || x2 > x + radius)
				return false;
			if(y2 < y - radius || y2 > y + radius)
				return false;
			return true;
	}

	public boolean isDragged() {
		return dragged;
	}

	public void setDragged(boolean dragged) {
		this.dragged = dragged;
	}
	
	public ImagePlus getImage(ImagePlus container, int size)
	{
		int radius = size/2;
		Rectangle r = new Rectangle(x - radius , y - radius, radius * 2, radius * 2);
		container.setRoi(r);
		ImageProcessor processor = container.getProcessor().crop();
		return new ImagePlus("", processor);
	}
	
	public ImageIcon getImageIcon(ImagePlus container, int size)
	{
		ImagePlus img = getImage(container, size);
		ImageIcon icon = new ImageIcon(img.getImage());
		return icon;
	}

	public static boolean boxContainedOnImage(int x, int y, int size, ImagePlus img) {
		if(img == null)
			return false;
		int width = img.getWidth();
		int height = img.getHeight();
		int radius = size/2;
		if(x - radius < 0)
			return false;
		if(x + radius > width)
			return false;
		if(y - radius < 0)
			return false;
		if(y + radius > height)
			return false;
		return true;
			
	}

	
}
