/*
 * GreatCircleDiffusionModel.java
 *
 * Copyright (C) 2002-2009 Alexei Drummond and Andrew Rambaut
 *
 * This file is part of BEAST.
 * See the NOTICE file distributed with this work for additional
 * information regarding copyright ownership and licensing.
 *
 * BEAST is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * BEAST is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with BEAST; if not, write to the
 * Free Software Foundation, Inc., 51 Franklin St, Fifth Floor,
 * Boston, MA  02110-1301  USA
 */

package beast.evolution.geo;

import beast.core.Description;
import beast.core.Input;
import beast.core.Input.Validate;
import beast.core.parameter.RealParameter;
import beast.evolution.substitutionmodel.ContinuousSubstitutionModel;

@Description("Diffusion model that assumes a normal diffusion process on a sphere, but without crossing a set of lines")
public class ConstrainedSphericalDiffusionModel extends ContinuousSubstitutionModel {
	public Input<RealParameter> precisionInput = new Input<RealParameter>("precision", "precision of diffusion process", Validate.REQUIRED);
	public Input<LineSet> tabuLineInput = new Input<LineSet>("tabulines", "set of lines that are not allowed to be crossed");

	RealParameter precision;

	double [][] latitudes;
	double [][] longitudes;
	
	
    @Override
    public void initAndValidate() throws Exception {
    	precision = precisionInput.get();
    	latitudes = tabuLineInput.get().getLatitudes();
    	longitudes = tabuLineInput.get().getLongitudes();
    	for (int i = 0; i < longitudes.length; i++) {
    		if (longitudes[i][0] < 0) longitudes[i][0] += 360;
    		if (longitudes[i][1] < 0) longitudes[i][1] += 360;
    	}
    	
    	super.initAndValidate();
    }

	@Override
    public double getLogLikelihood(double[] start, double[] stop, double time) {
		
		if (time <= 1e-20) {
			return -1e100;
		}
		if (start[0] == stop[0] && start[1] == stop[1]) {
			return -1e100;
		}
		double latitude1 = start[0];
		double longitude1 = start[1];
		double latitude2 = stop[0];
		double longitude2 = stop[1];
		if (longitude1 < 0) longitude1 += 360;
		if (longitude2 < 0) longitude2 += 360;
		
		// don't go through the atlantic
		if ((longitude1 < 330 && longitude1 > 210 && (longitude2 > 330 || longitude2 < 90)) ||
			(longitude2 < 330 && longitude2 > 210 && (longitude1 > 330 || longitude1 < 90))) {
			return -1e20 * (1 + 90 - latitude1 + 90 - latitude2);
		}
		// don't go through the pacific, unless over 45 latitude
		if ((latitude1 < 45 || latitude2 < 45) &&
			((longitude1 < 225 && longitude1 > 90 && longitude2 > 225 && longitude2 < 360) ||
			 (longitude2 < 225 && longitude2 > 90 && longitude1 > 225 && longitude1 < 360))) {
			return -1e20 * (1 + 90 - latitude1 + 90 - latitude2);
		}
		// don't go through the poles
		if (Math.abs(latitude1) > 70) {
			return -1e6 * (1+Math.abs(latitude1) - 70);
		}
		if (Math.abs(latitude2) > 70) {
			return -1e6 * (1+Math.abs(latitude2) - 70);
		}
		
		
		// check the line does not cross any of the tabu lines
//		for (int i = 0; i < latitudes.length; i++) {
//			if (intersects(latitude1, longitude1, latitude2, longitude2,
//					latitudes[i][0], longitudes[i][0], latitudes[i][1], longitudes[i][1])) {
//				return -1e100;
//			}
//		}
		
		
		// assumes start = {latitude, longitude}
		// assumes stop = {latitude, longitude}
		// and -90 < latitude < 90, -180 < longitude < 180
		
		double theta1 = (latitude1)*Math.PI/180.0;
		double phi1 = longitude1 * Math.PI/180;

		double theta2 = (latitude2)*Math.PI/180.0;
		double phi2 = longitude2 * Math.PI/180;
		
		double Deltalambda = phi2 - phi1;
		
		double angle = Math.acos(Math.sin(theta1)*Math.sin(theta2)+Math.cos(theta1) * Math.cos(theta2) * Math.cos(Deltalambda)); 

        double inverseVariance = precision.getValue(0) / time;
        double logP = -angle*angle * inverseVariance /2.0 + 0.5 * Math.log(angle * Math.sin(angle) * inverseVariance);
//		System.err.println(start[0] + " " + start[1] + " -> " + stop[0] + " " + stop[1] + " => " + logP);
        return logP;
        
    }

	private boolean intersects(double y1, double x1, double y2, double x2, 
			double y3, double x3, double y4, double x4) {
		double x=((x1*y2-y1*x2)*(x3-x4)-(x1-x2)*(x3*y4-y3*x4))/((x1-x2)*(y3-y4)-(y1-y2)*(x3-x4));
		double y=((x1*y2-y1*x2)*(y3-y4)-(y1-y2)*(x3*y4-y3*x4))/((x1-x2)*(y3-y4)-(y1-y2)*(x3-x4));
	    if (Double.isNaN(x)||Double.isNaN(y)) {
	        return false;
	    } else {
	        if (x1>=x2) {
	            if (!(x2<=x&&x<=x1)) {return false;}
	        } else {
	            if (!(x1<=x&&x<=x2)) {return false;}
	        }
	        if (y1>=y2) {
	            if (!(y2<=y&&y<=y1)) {return false;}
	        } else {
	            if (!(y1<=y&&y<=y2)) {return false;}
	        }
	        if (x3>=x4) {
	            if (!(x4<=x&&x<=x3)) {return false;}
	        } else {
	            if (!(x3<=x&&x<=x4)) {return false;}
	        }
	        if (y3>=y4) {
	            if (!(y4<=y&&y<=y3)) {return false;}
	        } else {
	            if (!(y3<=y&&y<=y4)) {return false;}
	        }
	    }
	    return true;	}

	public static void main(String[] args) {
		
		double [] start = new double[]{90, 0};
		double [] stop= new double[]{-90, 0};
		
		double latitude1 = start[0];
		double longitude1 = start[1];
		double theta1 = (latitude1)*Math.PI/180.0;
		if (longitude1 < 0) longitude1 += 360;
		double phi1 = longitude1 * Math.PI/180;

		double latitude2 = stop[0];
		double longitude2 = stop[1];
		double theta2 = (latitude2)*Math.PI/180.0;
		if (longitude2 < 0) longitude2 += 360;
		double phi2 = longitude2 * Math.PI/180;
		
		double Deltalambda = phi2 - phi1;
		
		double f1 = Math.sin(theta1)*Math.sin(theta2);
		double f2 = Math.cos(theta1) * Math.cos(theta2) * Math.cos(Deltalambda);
		
		double angle = Math.acos(f1+ f2);
		System.err.println("angle = " + angle * 180 / Math.PI);
		
	}

}
