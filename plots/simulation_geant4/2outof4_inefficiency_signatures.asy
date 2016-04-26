import root;
import pad_layout;

string topDir = "../../simulation_geant4/";

//string f = topDir + "analysis_full.root";
string f = topDir + "analysis.root";

string arms[] = { "L", "R" };

string pairs[] = { "T", "B" };

TH2_palette = Gradient(white, blue, heavygreen, yellow, red);

//----------------------------------------------------------------------------------------------------

for (int ai : arms.keys)
{
	for (int pi : pairs.keys)
	{
		string dir = arms[ai] + "/acs_" + pairs[pi] + "/";

		NewRow();

		NewPad(false);
		label(arms[ai] + ", " + pairs[pi]);

		NewPad("mean number of clusters $\bar N_c$");
		draw(RootGetObject(f, dir+"/h_N"), "vl", black, "near");
		draw(RootGetObject(f, dir+"/h_F"), "vl", red, "far");
		draw(RootGetObject(f, dir+"/h_NO"), "vl", blue, "opposite near");
		draw(RootGetObject(f, dir+"/h_FO"), "vl", heavygreen, "opposite far");
		limits((0, 0), (30, 1000), Crop);
		AttachLegend();

		NewPad("far: $\bar N_c$", "near: $\bar N_c$", axesAbove=true);
		draw(RootGetObject(f, dir+"/h2_N_F"));

		NewPad("far oppos.: $\bar N_c$", "near: $\bar N_c$", axesAbove=true);
		draw(RootGetObject(f, dir+"/h2_N_FO"));

		NewPad("far oppos.: $\bar N_c$", "far: $\bar N_c$", axesAbove=true);
		draw(RootGetObject(f, dir+"/h2_F_FO"));
		draw((0, 0)--(30, 30), black+1pt);
	}
}
