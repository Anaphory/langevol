package beast.evolution.geo;

import beast.core.BEASTObject;
import beast.core.Description;
import beast.core.Input;
import beast.core.Input.Validate;

@Description("Line, part of a LineSet")
public class Line extends BEASTObject {
	public Input<Double> xInput = new Input<Double>("x", "x coordinate of the line", Validate.REQUIRED);
	public Input<Double> yInput = new Input<Double>("y", "y coordinate of the line", Validate.REQUIRED);
	public Input<Double> wInput = new Input<Double>("w", "width of the line, x2 = x + w", Validate.REQUIRED);
	public Input<Double> hInput = new Input<Double>("h", "height of the line, y2 = y + w", Validate.REQUIRED);

	@Override
	public void initAndValidate() throws Exception {
	}

}
