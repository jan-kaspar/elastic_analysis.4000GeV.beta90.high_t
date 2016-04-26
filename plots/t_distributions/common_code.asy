
//----------------------------------------------------------------------------------------------------

real A_ref_8TeV = 530, B_ref_8TeV = 19.6;		// 8 TeV reference
real A_ref_7TeV = 506.51, B_ref_7TeV = 19.894;	// 7 TeV (trilogy) reference (fit 0.005 to 0.2 GeV^2)

//real A_ref = A_ref_8TeV, B_ref = B_ref_8TeV;
real A_ref = A_ref_7TeV, B_ref = B_ref_7TeV;

string MakeRefStr(string label="")
{
	string s = format("%.1f\,", A_ref) + format("\e^{-%.2f\,|t|}", B_ref);
	if (label != "")
		s += "\hbox{ ("+label+")}";
	return s;
}

string ref_str = MakeRefStr("");

//----------------------------------------------------------------------------------------------------

void DrawRelDiff(transform t = identity(), RootObject o, pen p, marker m = mCi+2pt+black, string label="")
{
	if (o.InheritsFrom("TH1"))
	{
		int N = o.iExec("GetNbinsX");
		for (int i = 1; i <= N; ++i)
		{
			real xc = o.rExec("GetBinCenter", i);
			real xw = o.rExec("GetBinWidth", i);
	
			real y = o.rExec("GetBinContent", i);
			real y_unc = o.rExec("GetBinError", i);
	
			real y_ref = A_ref * exp(-B_ref * xc);
	
			real y_rel = (y - y_ref) / y_ref;
			real y_rel_unc = y_unc / y_ref;
	
			draw(t * ((xc-xw/2, y_rel)--(xc+xw/2, y_rel)), p+0.1pt);	
			draw(t * ((xc, y_rel-y_rel_unc)--(xc, y_rel+y_rel_unc)), p+0.1pt);	
			draw(t * (xc, y_rel), mCi+0.001pt+p);
		}
	
		if (label != "")
			AddToLegend(label, mPl+4pt+p);
	}

	if (o.InheritsFrom("TGraph"))
	{
		guide g_u, g_b;
		
		int N = o.iExec("GetN");
		for (int i = 0; i <= N; ++i)
		{
			real xa[] = {0.};
			real ya[] = {0.};
			o.vExec("GetPoint", i, xa, ya);
			real x = xa[0];
			real y = ya[0];
	
			real y_unc = o.rExec("GetErrorY", i);
	
			real y_ref = A_ref * exp(-B_ref * x);
	
			real y_rel = (y - y_ref) / y_ref;
			real y_rel_unc = y_unc / y_ref;
	
			draw(t * ((x, y_rel-y_rel_unc)--(x, y_rel+y_rel_unc)), p);	
			draw(t * (x, y_rel), m);
			//g_u = g_u -- (x, y_rel+y_rel_unc);
			//g_b = g_b -- (x, y_rel-y_rel_unc);
		}

		//g_b = reverse(g_b);
		//filldraw(t*(g_u--g_b--cycle), p, nullpen);
		
		if (label != "")
			AddToLegend(label, mPl+5pt+p);
	}

	if (o.InheritsFrom("TF1"))
	{
		guide g;

		//real t_max = o.rExec("GetXmax");
		real t_max = 0.2;

		for (real t=0; t <= t_max; t += 0.001)
		{
			real y = o.rExec("Eval", t);
			real y_ref = A_ref * exp(-B_ref * t);
			real y_rel = (y - y_ref) / y_ref;
			g = g--(t, y_rel);
		}

		draw(g, p);
		
		if (label != "")
			AddToLegend(label, p);
	}
}

//----------------------------------------------------------------------------------------------------

void DrawRelDiffWithCorrection(transform t = identity(), RootObject o, RootObject o_c, real fc, pen p, string label="",
	bool arrow=false)
{
	if (o.InheritsFrom("TH1") && o_c.InheritsFrom("TH1"))
	{
		int N = o.iExec("GetNbinsX");
		for (int i = 1; i <= N; ++i)
		{
			real xc = o.rExec("GetBinCenter", i);
			real xw = o.rExec("GetBinWidth", i);
	
			real y = o.rExec("GetBinContent", i);
			real y_unc = o.rExec("GetBinError", i);
			
			real xc_corr = o_c.rExec("GetBinCenter", i);
			real corr = o_c.rExec("GetBinContent", i);

			if (abs(xc_corr - xc) > 0.001)
				write("ERROR: binning mismatch");

			// scale correction
			corr = 1. + fc * (corr - 1.);
	
			// apply correction
			real y_corr = y * corr;
			real y_unc_corr = y_unc * corr;
	
			real y_ref = A_ref * exp(-B_ref * xc);
	
			real y_rel = (y - y_ref) / y_ref;
			real y_rel_unc = y_unc / y_ref;

			real y_corr_rel = (y_corr - y_ref) / y_ref;
			real y_unc_corr_rel = y_unc_corr / y_ref;
	
			if (arrow)
			{
				draw(t * ((xc, y_rel)--(xc, y_corr_rel)), p, EndArrow(5));		
			} else {
				draw(t * ((xc-xw/2, y_corr_rel)--(xc+xw/2, y_corr_rel)), p+0.1pt);	
				draw(t * ((xc, y_corr_rel-y_unc_corr_rel)--(xc, y_corr_rel+y_unc_corr_rel)), p+0.1pt);	
				//draw(t * (xc, y_rel_corr), mCi+0.001pt+p);
			}
		}

		if (label != "")
			AddToLegend(label, mPl+4pt+p);
	}

	if (o.InheritsFrom("TH1") && o_c.InheritsFrom("TF1"))
	{
		int N = o.iExec("GetNbinsX");
		for (int i = 1; i <= N; ++i)
		{
			real xc = o.rExec("GetBinCenter", i);
			real xw = o.rExec("GetBinWidth", i);
	
			real y = o.rExec("GetBinContent", i);
			real y_unc = o.rExec("GetBinError", i);
			
			real corr = o_c.rExec("Eval", xc);

			// scale correction
			corr = 1. + fc * (corr - 1.);
	
			// apply correction
			real y_corr = y * corr;
			real y_unc_corr = y_unc * corr;
	
			real y_ref = A_ref * exp(-B_ref * xc);
	
			real y_rel = (y - y_ref) / y_ref;
			real y_rel_unc = y_unc / y_ref;

			real y_corr_rel = (y_corr - y_ref) / y_ref;
			real y_unc_corr_rel = y_unc_corr / y_ref;
	
			if (arrow)
			{
				draw(t * ((xc, y_rel)--(xc, y_corr_rel)), p, EndArrow(5));		
			} else {
				draw(t * ((xc-xw/2, y_corr_rel)--(xc+xw/2, y_corr_rel)), p+0.1pt);	
				draw(t * ((xc, y_corr_rel-y_unc_corr_rel)--(xc, y_corr_rel+y_unc_corr_rel)), p+0.1pt);	
				//draw(t * (xc, y_rel_corr), mCi+0.001pt+p);
			}
		}

		if (label != "")
			AddToLegend(label, mPl+4pt+p);
	}
}

//----------------------------------------------------------------------------------------------------

void DrawRelDiffWithCorrectionArray(transform t = identity(), RootObject o, real corr_a[], pen p, string label="",
	bool arrow=false)
{
	if (o.InheritsFrom("TH1"))
	{
		int N = o.iExec("GetNbinsX");
		for (int i = 1; i <= N; ++i)
		{
			real xc = o.rExec("GetBinCenter", i);
			real xw = o.rExec("GetBinWidth", i);
	
			real y = o.rExec("GetBinContent", i);
			real y_unc = o.rExec("GetBinError", i);
			
			real corr = corr_a[i];

			// apply correction
			real y_corr = y * corr;
			real y_unc_corr = y_unc * corr;
	
			real y_ref = A_ref * exp(-B_ref * xc);
	
			real y_rel = (y - y_ref) / y_ref;
			real y_rel_unc = y_unc / y_ref;

			real y_corr_rel = (y_corr - y_ref) / y_ref;
			real y_unc_corr_rel = y_unc_corr / y_ref;
	
			if (arrow)
			{
				draw(t * ((xc, y_rel)--(xc, y_corr_rel)), p+0.1pt, EndArrow(5));		
			} else {
				draw(t * ((xc-xw/2, y_corr_rel)--(xc+xw/2, y_corr_rel)), p+0.1pt);	
				draw(t * ((xc, y_corr_rel-y_unc_corr_rel)--(xc, y_corr_rel+y_unc_corr_rel)), p+0.1pt);	
				//draw(t * (xc, y_rel_corr), mCi+0.001pt+p);
			}
		}
	}

	AddToLegend(label, mPl+4pt+p);
}

//----------------------------------------------------------------------------------------------------

void DrawRelDiffWithBand(transform t = identity(), RootObject o, RootObject o_u, real unc_scale = 1., pen p, string label)
{
	if (o.InheritsFrom("TH1"))
	{
		guide g_u, g_b;
		guide g_u_full, g_b_full;

		int N = o.iExec("GetNbinsX");
		for (int i = 1; i <= N; ++i)
		{
			real xc = o.rExec("GetBinCenter", i);
			real xw = o.rExec("GetBinWidth", i);
	
			real y = o.rExec("GetBinContent", i);
			real y_unc = o.rExec("GetBinError", i);
	
			real y_ref = 530 * exp(-19.6 * xc);
	
			real y_rel = (y - y_ref) / y_ref;
			real y_rel_unc = y_unc / y_ref;
	
			draw(t * ((xc-xw/2, y_rel)--(xc+xw/2, y_rel)), p+0.1pt);	
			draw(t * ((xc, y_rel-y_rel_unc)--(xc, y_rel+y_rel_unc)), p+0.1pt);	
			draw(t * (xc, y_rel), mCi+1pt+p);

			real y_rel_sys_unc = unc_scale * (o_u.rExec("GetBinContent", i) - 0.) * y / y_ref;

			g_u = g_u -- (xc, y_rel+y_rel_sys_unc);
			g_b = g_b -- (xc, y_rel-y_rel_sys_unc);

			y_rel_sys_unc = sqrt(0.04*0.04 + y_rel_sys_unc^2);

			g_u_full = g_u_full -- (xc, y_rel+y_rel_sys_unc);
			g_b_full = g_b_full -- (xc, y_rel-y_rel_sys_unc);
		}

		//g_b_full = reverse(g_b_full);
		//filldraw(t*(g_u_full--g_b_full--cycle), heavygreen+opacity(0.3), nullpen);

		g_b = reverse(g_b);
		filldraw(t*(g_u--g_b--cycle), p+opacity(0.3), nullpen);
	}

	AddToLegend(label, p);
}


//----------------------------------------------------------------------------------------------------

void DrawRelDiffBand(transform t = identity(), RootObject o, RootObject o_u,
	real x_max = +inf,
	real unc_rel_const = 0., real unc_scale = 1.,
	pen p, string label)
{
	if (o.InheritsFrom("TH1"))
	{
		guide g_u, g_b;
		guide g_u_full, g_b_full;

		int N = o.iExec("GetNbinsX");
		for (int i = 1; i <= N; ++i)
		{
			real xc = o.rExec("GetBinCenter", i);
			real xw = o.rExec("GetBinWidth", i);
	
			real y = o.rExec("GetBinContent", i);
			real y_stat_unc = o.rExec("GetBinError", i);

			if (fabs(y) < 1e-3)
				continue;

			if (xc > x_max)
				continue;
	
			real y_ref = A_ref * exp(-B_ref * xc);
	
			real y_rel = (y - y_ref) / y_ref;
			real y_rel_stat_unc = y_stat_unc / y_ref;
	
			real y_rel_sys_unc = o_u.rExec("GetBinContent", i);
			y_rel_sys_unc = unc_scale * sqrt(y_rel_sys_unc^2 + unc_rel_const^2);

			y_rel_sys_unc = y_rel_sys_unc * y / y_ref;

			real y_band_cen = 5.43464e+02 * exp( -2.07226e+01*xc + 1.03308e+01*xc^2 -2.33729e+01*xc^3 );
			real y_rel_band_cen = (y_band_cen - y_ref) / y_ref;

			g_u = g_u -- (xc, y_rel_band_cen + y_rel_sys_unc);
			g_b = g_b -- (xc, y_rel_band_cen - y_rel_sys_unc);
		}

		g_b = reverse(g_b);
		filldraw(t*(g_u--g_b--cycle), p+opacity(1), nullpen);
	
		if (label != "")
			AddToLegend(label, mSq+6pt+p);
	}
	
	if (o.InheritsFrom("TF1"))
	{
		guide g_u, g_b;
		guide g_u_full, g_b_full;

		int N = o_u.iExec("GetNbinsX");
		for (int bi = 1; bi <= N; ++bi)
		{
			real xc = o_u.rExec("GetBinCenter", bi);
			real y_rel_unc = o_u.rExec("GetBinContent", bi);

			if (y_rel_unc == 0)
				continue;

			real y_cen = o.rExec("Eval", xc);

			if (xc > x_max)
				continue;
			
			real y_unc = y_rel_unc * y_cen;
	
			real y_ref = A_ref * exp(-B_ref * xc);
	
			real y_cen_rel = (y_cen - y_ref) / y_ref;
			real y_unc_rel = y_unc / y_ref;

			g_u = g_u -- (xc, y_cen_rel + y_unc_rel);
			g_b = g_b -- (xc, y_cen_rel - y_unc_rel);
		}

		g_b = reverse(g_b);
		filldraw(t*(g_u--g_b--cycle), p, nullpen);

		if (label != "")
			AddToLegend(label, mSq+6pt+p);
	}
}

//----------------------------------------------------------------------------------------------------

void DrawExpFit(real A, real B, pen p, string label="")
{
	guide g;

	for (real t = 0; t <= 0.3; t += 0.003)
	{
		real y = A * exp(-B * t);
		real y_ref = A_ref * exp(-B_ref * t);

		real y_rel = (y - y_ref) / y_ref;

		g = g -- (t, y_rel);
	}

	draw(g, p, label);
}
