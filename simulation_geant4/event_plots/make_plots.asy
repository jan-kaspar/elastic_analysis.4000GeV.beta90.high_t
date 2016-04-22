import pad_layout;

//----------------------------------------------------------------------------------------------------

struct Particle
{
	int code;
	string name;
	pen p;
}

Particle particles[];

void AddParticle(int _c, string _n, pen _p)
{
	Particle p;
	p.code = _c;
	p.name = _n;
	p.p = _p;
	particles.push(p);
}

AddParticle(2212, "p", red);
AddParticle(22, "$\ga$", heavygreen);
AddParticle(+11, "$\rm e^-$", blue);
AddParticle(-11, "$\rm e^+$", cyan);
AddParticle(+211, "$\rm\pi^+$", olive);
AddParticle(-211, "$\rm\pi^-$", magenta);

pen GetParticlePen(int code)
{
	for (int pi : particles.keys)
		if (particles[pi].code == code)
			return particles[pi].p;

	return black;
}

//----------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------

struct Zone
{
	string name;
	real z0;
	real edge_top, edge_bot;
}

Zone zones[];

void AddZone(string _n, real _z0, real _et, real _eb)
{
	Zone z;
	z.name = _n;
	z.z0 = _z0;
	z.edge_top = _et;
	z.edge_bot = _eb;

	zones.push(z);
}

AddZone("Left Far", -220.000e3, 8.0, -8.2);
AddZone("Left Near", -214.628e3, 7.2, -7.4);
AddZone("Right Near", 214.628e3, 7.2, -7.4);
AddZone("Right Far", 220.000e3, 8.0, -8.2);

//----------------------------------------------------------------------------------------------------

pad zone_pads[];

real z_array[] = { -20.3, -15.7, -11.3, -6.7, -2.3, 2.3, 6.7, 11.3, 15.7, 20.3 };

string event_desc;

void BeginEvent(int id, string _ed)
{
	write(">> event ", id);

	NewRow();

	event_desc = _ed;

	zone_pads.delete();
	for (int zi : zones.keys)
	{
		pad p =	NewPad(format("$z - (%.0f)\ung{mm}$", zones[zi].z0), "$y\ung{mm}$");

		for (int i = 0; i < 10; ++i)
		{
			real z = z_array[i];

			real gap = 0.35, len = 39.88;

			draw((z, zones[zi].edge_top + gap)--(z, zones[zi].edge_top + gap + len), lightgray+1pt);
			draw((z, zones[zi].edge_bot - gap)--(z, zones[zi].edge_bot - gap - len), lightgray+1pt);
		}

		draw((-30, zones[zi].edge_top)--(+30, zones[zi].edge_top), lightgray+1pt);
		draw((-30, zones[zi].edge_bot)--(+30, zones[zi].edge_bot), lightgray+1pt);

		zone_pads.push(p);
	}
}

//----------------------------------------------------------------------------------------------------

int GetZone(real z)
{
	for (int zi : zones.keys)
	{
		real th = 100;	// mm
		if (fabs(z - zones[zi].z0) < th)
			return zi;
	}

	return -1;
}

//----------------------------------------------------------------------------------------------------

int GetZone(int rpId)
{
	if (rpId == 24 || rpId == 25) return 0;
	if (rpId == 20 || rpId == 21) return 1;
	if (rpId == 120 || rpId == 121) return 2;
	if (rpId == 124 || rpId == 125) return 3;

	return -1;
}

//----------------------------------------------------------------------------------------------------

void Hit(real x, real y, real z, real energyLoss, int track, int particle)
{
	pen p = GetParticlePen(particle);
	
	// determine zone
	int zi = GetZone(z);

	if (zi < 0)
		return;

	// draw point
	SetPad(zone_pads[zi]);
	dot((z - zones[zi].z0, y), p);	
}

//----------------------------------------------------------------------------------------------------

void TrackFit(int rpId, real x, real y, real z, real th_x, real th_y, bool full_sim)
{
	int zi = GetZone(z);
	if (zi < 0)
		return;

	real z0 = zones[zi].z0;

	real de_z = 40;
	pair l = (z - z0, y) - de_z * (1., th_y);
	pair r = (z - z0, y) + de_z * (1., th_y);
	pen p = (full_sim) ? black : black + dashed;

	SetPad(zone_pads[zi]);
	draw(l--r, p);
}

//----------------------------------------------------------------------------------------------------

void RPClusters(int rpId, bool full_sim ... int occupancy[])
{
	int zi = GetZone(rpId);
	if (zi < 0)
		return;

	real y_sign = ((rpId % 2) == 0) ? +1 : -1;
	real z_sign = (quotient(rpId, 100) > 0) ? +1 : -1;

	real y_abs = 40;

	SetPad(zone_pads[zi]);

	for (int i : z_array.keys)
	{
		label(format("%u", occupancy[i]), (z_sign * z_array[i], y_sign*y_abs));
	}
}

//----------------------------------------------------------------------------------------------------

void EndEvent(int id)
{
	for (int zi : zones.keys)
	{
		SetPad(zone_pads[zi]);
		limits((-40, -50), (+40, +50), Crop);
		AttachLegend(zones[zi].name, E, E);
	}

	NewPad(false);
	AddToLegend("<{\it Geant4}:");
	for (int pi : particles.keys)
	{
		AddToLegend(particles[pi].name, mCi+2pt+particles[pi].p);
	}
	AddToLegend("other", mCi+2pt+black);

	AddToLegend("<{\it reco tracks}:");
	AddToLegend("full reco", black);
	AddToLegend("ideal reco", black+dashed);

	AttachLegend(format("event %u, ", id) + replace(replace(event_desc, "&", "\&"), "_", "\_"), E, E);

	//GShipout(format("event_%06u", id), hSkip=1mm);	
}

//----------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------

include "dump_data.asy";
