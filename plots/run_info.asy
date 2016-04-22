//----------------------------------------------------------------------------------------------------

transform swToHours = scale(1/3600, 1);

//----------------------------------------------------------------------------------------------------

int runs[];
real ts_from[], ts_to[];
pen colors[];

void AddRun(int r, real f, real t, pen p = yellow)
{
	runs.push(r);
	ts_from.push(f);
	ts_to.push(t);
	colors.push(p+opacity(0.3));
}

AddRun(8369, 71048, 83218);
AddRun(8371, 83777, 95420);
AddRun(8372, 95804, 112410);

//----------------------------------------------------------------------------------------------------

void DrawRunBands(transform t = identity(), real y_min=0, real y_max=0, bool details=true)
{
	for (int i : runs.keys)
	{
		real x_min = ts_from[i]/3600, x_max = ts_to[i]/3600;

		pen p = (details) ? colors[i]+opacity(0.3) : yellow+opacity(0.3);
		filldraw(t * ((x_min, y_min)--(x_max, y_min)--(x_max, y_max)--(x_min, y_max)--cycle), p, nullpen);

		if (details)
			label(format("{\SmallerFonts %u}", runs[i]), t * ((x_min + x_max)/2, y_max), S);
	}
}

//----------------------------------------------------------------------------------------------------

void DrawRunBoundaries()
{
	for (int i : runs.keys)
	{
		yaxis(XEquals(ts_from[i]/3600, false), dashed);
		yaxis(XEquals(ts_to[i]/3600, false), dashed);
	}
}
