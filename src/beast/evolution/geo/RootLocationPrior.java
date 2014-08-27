package beast.evolution.geo;

import java.util.List;
import java.util.Random;

import beast.core.Description;
import beast.core.Distribution;
import beast.core.Input;
import beast.core.Input.Validate;
import beast.core.State;
import beast.core.parameter.RealParameter;

@Description("Root prior putting the location in Africa")
public class RootLocationPrior extends Distribution {

	public Input<RealParameter> locationInput = new Input<RealParameter>("trait", "location trait where last two elements represent latitude, longitude", Validate.REQUIRED);
	
	RealParameter location;
	
	@Override
	public void initAndValidate() throws Exception {
		location = locationInput.get();
	}
	
	@Override
	public double calculateLogP() throws Exception {
		logP = 0;
		int i = location.getDimension() - 2;
		double latitude = location.getArrayValue(i);
		double longitude = location.getArrayValue(i+1);
		if (latitude > 22) {
			logP += latitude * -1e10;
		}
		if (longitude > 38) {
			logP += longitude * -1e10;
		}
		return logP;
	}
	
	
	
	@Override
	public List<String> getArguments() {return null;}
	@Override
	public List<String> getConditions() {return null;}
	@Override
	public void sample(State state, Random random) {}

}
