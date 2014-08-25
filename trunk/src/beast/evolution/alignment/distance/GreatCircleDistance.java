package beast.evolution.alignment.distance;

import beast.core.BEASTObject;
import beast.core.Description;
import beast.core.Input;
import beast.core.Input.Validate;
import beast.evolution.tree.TreeInterface;
import beast.evolution.tree.TreeTraitMap;

@Description("Distance between points on a globe")
public class GreatCircleDistance extends BEASTObject implements Distance {
	public Input<TreeTraitMap> traitInput = new Input<TreeTraitMap>("trait", "trait specifying latitude/longitude locations", Validate.REQUIRED);
			
	TreeTraitMap trait;
	TreeInterface tree;

	final static double EARTHRADIUS = 6371; // mean radius, according to http://en.wikipedia.org/wiki/Earth_radius
	
	@Override
	public void initAndValidate() throws Exception {
		trait = traitInput.get();
		tree = trait.treeInput.get();
	}

	@Override
	public double pairwiseDistance(int taxon1, int taxon2) {
		double [] loc1 = trait.getTrait(tree, tree.getNode(taxon1));
		double [] loc2 = trait.getTrait(tree, tree.getNode(taxon2));
		
		double latitude1 = loc1[0];
		double longitude1 = loc1[1];
		double theta1 = (latitude1)*Math.PI/180.0;
		if (longitude1 < 0) longitude1 += 360;
		double phi1 = longitude1 * Math.PI/180;

		double latitude2 = loc2[0];
		double longitude2 = loc2[1];
		double theta2 = (latitude2)*Math.PI/180.0;
		if (longitude2 < 0) longitude2 += 360;
		double phi2 = longitude2 * Math.PI/180;
		
		double Deltalambda = phi2 - phi1;
		
		double angle = Math.acos(Math.sin(theta1)*Math.sin(theta2)+Math.cos(theta1) * Math.cos(theta2) * Math.cos(Deltalambda)); 

		return angle * EARTHRADIUS;
	}

}
