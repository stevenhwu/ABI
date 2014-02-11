package dr.app.beauti.components.continuous;

import dr.app.beauti.components.ComponentFactory;
import dr.app.beauti.generator.ComponentGenerator;
import dr.app.beauti.options.BeautiOptions;
import dr.app.beauti.options.ComponentOptions;

/**
 * @author Andrew Rambaut
 * @version $Id: ContinuousComponentFactory.java 5861 2013-10-03 13:38:22Z rambaut $
 */

public class ContinuousComponentFactory implements ComponentFactory {

	public static ComponentFactory INSTANCE = new ContinuousComponentFactory();

	public ComponentGenerator createGenerator(BeautiOptions beautiOptions) {
		return new ContinuousComponentGenerator(beautiOptions);
	}

	public ComponentOptions createOptions(BeautiOptions beautiOptions) {
        return new ContinuousComponentOptions(beautiOptions);
	}

}