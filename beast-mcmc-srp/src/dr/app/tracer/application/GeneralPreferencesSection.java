package dr.app.tracer.application;

import jam.preferences.PreferencesSection;
import jam.util.IconUtils;

import javax.swing.*;
import java.util.prefs.Preferences;

/**
 * @author Andrew Rambaut
 * @version $Id: GeneralPreferencesSection.java,v 1.1 2006/09/09 15:23:33 rambaut Exp $
 */
public class GeneralPreferencesSection implements PreferencesSection {
	Icon projectToolIcon = IconUtils.getIcon(TracerApp.class, "images/prefsGeneral.png");

	public String getTitle() {
		return "General";
	}

	public Icon getIcon() {
		return projectToolIcon;
	}

	public JPanel getPanel() {
		JPanel panel = new JPanel();
		panel.add(generalCheck);
		return panel;
	}

	public void retrievePreferences() {
		Preferences prefs = Preferences.userNodeForPackage(TracerApp.class);
		generalCheck.setSelected(prefs.getBoolean("general_check", true));
	}

	public void storePreferences() {
		Preferences prefs = Preferences.userNodeForPackage(TracerApp.class);
		prefs.putBoolean("general_check", generalCheck.isSelected());
	}

	JCheckBox generalCheck = new JCheckBox("The preferences window is not implemented yet.");
}
