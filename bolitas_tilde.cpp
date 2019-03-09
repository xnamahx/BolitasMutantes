#include "c74_msp.h"
#include <iostream>

#include "marbles/random/random_generator.h"
#include "marbles/random/random_stream.h"
#include "marbles/random/t_generator.h"
#include "marbles/random/x_y_generator.h"
#include "marbles/note_filter.h"

using namespace c74::max;

static t_class* this_class = nullptr;

inline double constrain(double v, double vMin, double vMax) {
	return std::max<double>(vMin, std::min<double>(vMax, v));
}

static const int BLOCK_SIZE = 5;

static const marbles::Scale preset_scales[6] = {
	// C major
	{
		1.0f,
		12,
		{
			{ 0.0000f, 255 },  // C
			{ 0.0833f, 16 },   // C#
			{ 0.1667f, 96 },   // D
			{ 0.2500f, 24 },   // D#
			{ 0.3333f, 128 },  // E
			{ 0.4167f, 64 },   // F
			{ 0.5000f, 8 },    // F#
			{ 0.5833f, 192 },  // G
			{ 0.6667f, 16 },   // G#
			{ 0.7500f, 96 },   // A
			{ 0.8333f, 24 },   // A#
			{ 0.9167f, 128 },  // B
		}
	},

	// C minor
	{
		1.0f,
		12,
		{
			{ 0.0000f, 255 },  // C
			{ 0.0833f, 16 },   // C#
			{ 0.1667f, 96 },   // D
			{ 0.2500f, 128 },  // Eb
			{ 0.3333f, 8 },    // E
			{ 0.4167f, 64 },   // F
			{ 0.5000f, 4 },    // F#
			{ 0.5833f, 192 },  // G
			{ 0.6667f, 16 },   // G#
			{ 0.7500f, 96 },   // A
			{ 0.8333f, 128 },  // Bb
			{ 0.9167f, 16 },   // B
		}
	},

	// Pentatonic
	{
		1.0f,
		12,
		{
			{ 0.0000f, 255 },  // C
			{ 0.0833f, 4 },    // C#
			{ 0.1667f, 96 },   // D
			{ 0.2500f, 4 },    // Eb
			{ 0.3333f, 4 },    // E
			{ 0.4167f, 140 },  // F
			{ 0.5000f, 4 },    // F#
			{ 0.5833f, 192 },  // G
			{ 0.6667f, 4 },    // G#
			{ 0.7500f, 96 },   // A
			{ 0.8333f, 4 },    // Bb
			{ 0.9167f, 4 },    // B
		}
	},

	// Pelog
	{
		1.0f,
		7,
		{
			{ 0.0000f, 255 },  // C
			{ 0.1275f, 128 },  // Db+
			{ 0.2625f, 32 },  // Eb-
			{ 0.4600f, 8 },    // F#-
			{ 0.5883f, 192 },  // G
			{ 0.7067f, 64 },  // Ab
			{ 0.8817f, 16 },    // Bb+
		}
	},

	// Raag Bhairav That
	{
		1.0f,
		12,
		{
			{ 0.0000f, 255 }, // ** Sa
			{ 0.0752f, 128 }, // ** Komal Re
			{ 0.1699f, 4 },   //    Re
			{ 0.2630f, 4 },   //    Komal Ga
			{ 0.3219f, 128 }, // ** Ga
			{ 0.4150f, 64 },  // ** Ma
			{ 0.4918f, 4 },   //    Tivre Ma
			{ 0.5850f, 192 }, // ** Pa
			{ 0.6601f, 64 },  // ** Komal Dha
			{ 0.7549f, 4 },   //    Dha
			{ 0.8479f, 4 },   //    Komal Ni
			{ 0.9069f, 64 },  // ** Ni
		}
	},

	// Raag Shri
	{
		1.0f,
		12,
		{
			{ 0.0000f, 255 }, // ** Sa
			{ 0.0752f, 4 },   //    Komal Re
			{ 0.1699f, 128 }, // ** Re
			{ 0.2630f, 64 },  // ** Komal Ga
			{ 0.3219f, 4 },   //    Ga
			{ 0.4150f, 128 }, // ** Ma
			{ 0.4918f, 4 },   //    Tivre Ma
			{ 0.5850f, 192 }, // ** Pa
			{ 0.6601f, 4 },   //    Komal Dha
			{ 0.7549f, 64 },  // ** Dha
			{ 0.8479f, 128 }, // ** Komal Ni
			{ 0.9069f, 4 },   //    Ni
		}
	},
};

void* to1;
void* to2;
void* to3;
void* yo;
void* xo1;
void* xo2;
void* xo3;

marbles::RandomGenerator random_generator;
marbles::RandomStream random_stream;
marbles::TGenerator t_generator;
marbles::XYGenerator xy_generator;
marbles::NoteFilter note_filter;

// State
/*BooleanTrigger tDejaVuTrigger;
BooleanTrigger xDejaVuTrigger;
BooleanTrigger tModeTrigger;
BooleanTrigger xModeTrigger;
BooleanTrigger tRangeTrigger;
BooleanTrigger xRangeTrigger;
BooleanTrigger externalTrigger;*/
bool t_deja_vu;
bool x_deja_vu;
int t_mode;
int x_mode;
int t_range;
int x_range;
bool external;
int x_scale;
int y_divider_index;
int x_clock_source_internal;


float f_pulse_width_mean;
float f_pulse_width_std;
double f_t_bias;
double f_x_bias;
double f_deja_vu;
double f_deja_vu_length;
double deja_vu_length;
bool f_t_clock_input_patched;
bool f_t_clock_input;
bool f_x_clock_input_patched;
bool f_x_clock_input;

// Buffers
stmlib::GateFlags t_clocks[BLOCK_SIZE] = {};
stmlib::GateFlags last_t_clock = 0;
stmlib::GateFlags xy_clocks[BLOCK_SIZE] = {};
stmlib::GateFlags last_xy_clock = 0;
float ramp_master[BLOCK_SIZE] = {};
float ramp_external[BLOCK_SIZE] = {};
float ramp_slave[2][BLOCK_SIZE] = {};
bool gates[BLOCK_SIZE * 2] = {};
float voltages[BLOCK_SIZE * 4] = {};
int blockIndex = 0;

double f_t_rate;
double f_x_spread;

double f_t_jitter;
double f_x_steps;

int sampleRate;

static const int loop_length[] = {
	1, 1, 1, 2, 2,
	2, 2, 2, 3, 3,
	3, 3, 4, 4, 4,
	4, 4, 5, 5, 6,
	6, 6, 7, 7, 8,
	8, 8, 10, 10, 12,
	12, 12, 14, 14, 16,
	16
};

struct t_bolitas {
	t_object  x_obj;
};

void onReset() {
	t_deja_vu = false;
	x_deja_vu = false;
	t_mode = 0;
	x_mode = 0;
	t_range = 1;
	x_range = 2;
	external = false;
	f_t_clock_input_patched = false;
	f_x_clock_input_patched = false;
	f_t_clock_input = false;
	f_x_clock_input = false;
	x_scale = 0;
	y_divider_index = 8;
	x_clock_source_internal = 0;

	f_t_bias = 0.0f;
	f_x_bias = 0.0f;
	f_deja_vu = 0.0f;
	f_deja_vu_length = 0.0f;
	float deja_vu_length_index = f_deja_vu_length * (sizeof(loop_length) - 1);
	deja_vu_length = loop_length[(int) roundf(deja_vu_length_index)];

	f_t_jitter = 0.0f;
	f_pulse_width_std = 0.0f;
	f_pulse_width_mean = 0.0f;
}

void onSampleRateChange() {
	float sr = sys_getsr();
	if(sampleRate!=sr){
		sampleRate = sr;
		t_generator.Init(&random_stream, sampleRate);
		xy_generator.Init(&random_stream, sampleRate);

		// Set scales
		for (int i = 0; i < 6; i++) {
			xy_generator.LoadScale(i, preset_scales[i]);
		}
	}
}

void stepBlock() {
	// Ramps
	marbles::Ramps ramps;
	ramps.master = ramp_master;
	ramps.external = ramp_external;
	ramps.slave[0] = ramp_slave[0];
	ramps.slave[1] = ramp_slave[1];

	// Set up TGenerator
	bool t_external_clock = f_t_clock_input_patched;

	t_generator.set_model((marbles::TGeneratorModel) t_mode);
	t_generator.set_range((marbles::TGeneratorRange) t_range);
	double t_rate = 60.f * f_t_rate;
	t_generator.set_rate(t_rate);
	double t_bias = constrain(f_t_bias, 0.f, 1.f);
	t_generator.set_bias(t_bias);
	double t_jitter = constrain(f_t_jitter, 0.f, 1.f);
	t_generator.set_jitter(t_jitter);
	t_generator.set_deja_vu(t_deja_vu ? f_deja_vu : 0.f);
	t_generator.set_length(deja_vu_length);
	t_generator.set_pulse_width_mean(f_pulse_width_mean);
	t_generator.set_pulse_width_std(f_pulse_width_std);

	t_generator.Process(t_external_clock, t_clocks, ramps, gates, BLOCK_SIZE);

	// Set up XYGenerator
	marbles::ClockSource x_clock_source = (marbles::ClockSource) x_clock_source_internal;
	if (f_x_clock_input_patched)
		x_clock_source = marbles::CLOCK_SOURCE_EXTERNAL;

	marbles::GroupSettings x;
	x.control_mode = (marbles::ControlMode) x_mode;
	x.voltage_range = (marbles::VoltageRange) x_range;
	// TODO Fix the scaling
	double note_cv = 0.5f * (f_x_spread);
	double u = note_filter.Process(0.5f * (note_cv + 1.f));
	x.register_mode = external;
	x.register_value = u;

	double x_spread = constrain(f_x_spread, 0.f, 1.f);
	x.spread = x_spread;
	double x_bias = constrain(f_x_bias, 0.f, 1.f);
	x.bias = x_bias;
	double x_steps = constrain(f_x_steps, 0.f, 1.f);
	x.steps = x_steps;
	x.deja_vu = x_deja_vu ? f_deja_vu : 0.f;
	x.length = deja_vu_length;
	x.ratio.p = 1;
	x.ratio.q = 1;
	x.scale_index = x_scale;

	marbles::GroupSettings y;
	y.control_mode = marbles::CONTROL_MODE_IDENTICAL;
	// TODO
	y.voltage_range = (marbles::VoltageRange) x_range;
	y.register_mode = false;
	y.register_value = 0.0f;
	// TODO
	y.spread = x_spread;
	y.bias = x_bias;
	y.steps = x_steps;
	y.deja_vu = 0.0f;
	y.length = 1;
	static const marbles::Ratio y_divider_ratios[] = {
		{ 1, 64 },
		{ 1, 48 },
		{ 1, 32 },
		{ 1, 24 },
		{ 1, 16 },
		{ 1, 12 },
		{ 1, 8 },
		{ 1, 6 },
		{ 1, 4 },
		{ 1, 3 },
		{ 1, 2 },
		{ 1, 1 },
	};
	y.ratio = y_divider_ratios[y_divider_index];
	y.scale_index = x_scale;

	xy_generator.Process(x_clock_source, x, y, xy_clocks, ramps, voltages, BLOCK_SIZE);
}


void bolitas_run(t_bolitas *x, double f) {
		// Clocks
		last_t_clock = stmlib::ExtractGateFlags(last_t_clock, f_t_clock_input);
		t_clocks[blockIndex] = last_t_clock;

		last_xy_clock = stmlib::ExtractGateFlags(last_xy_clock, f_x_clock_input);
		xy_clocks[blockIndex] = last_xy_clock;

		// Process block
		if (++blockIndex >= BLOCK_SIZE) {
			//onSampleRateChange();
			blockIndex = 0;
			stepBlock();
		}

    	outlet_float(to1,gates[blockIndex*2 + 0] ? 10.f : 0.f);
		outlet_float(to2,(ramp_master[blockIndex] < 0.5f) ? 10.f : 0.f);
		outlet_float(to3,gates[blockIndex*2 + 1] ? 10.f : 0.f);
		outlet_float(xo1,voltages[blockIndex*4 + 0]);
		outlet_float(xo2,voltages[blockIndex*4 + 1]);
		outlet_float(xo3,voltages[blockIndex*4 + 2]);
		outlet_float(yo,voltages[blockIndex*4 + 3]);
}



void bolitas_perform64(t_bolitas* self, t_object* dsp64, double** ins, long numins, double** outs, long numouts, long sampleframes, long flags, void* userparam) {
    double    *t1 = outs[0];   // first outlet
    double    *t2 = outs[1];   // first outlet
    double    *t3 = outs[2];   // first outlet
    double    *y = outs[3];   // first outlet
    double    *x1 = outs[4];   // first outlet
    double    *x2 = outs[5];   // first outlet
    double    *x3 = outs[6];   // first outlet

    if (numouts>0)
    {
		// Clocks
		last_t_clock = stmlib::ExtractGateFlags(last_t_clock, f_t_clock_input);
		t_clocks[blockIndex] = last_t_clock;

		last_xy_clock = stmlib::ExtractGateFlags(last_xy_clock, f_x_clock_input);
		xy_clocks[blockIndex] = last_xy_clock;

		// Process block
		if (++blockIndex >= BLOCK_SIZE) {
			//onSampleRateChange();
			blockIndex = 0;
			stepBlock();
		}

    	*t1 = gates[blockIndex*2 + 0] ? 10.f : 0.f;
		*t2 = (ramp_master[blockIndex] < 0.5f) ? 10.f : 0.f;
		*t3 = gates[blockIndex*2 + 1] ? 10.f : 0.f;
		*x1 = voltages[blockIndex*4 + 0];
		*x2 = voltages[blockIndex*4 + 1];
		*x3 = voltages[blockIndex*4 + 2];
		*y = voltages[blockIndex*4 + 3];
	}
}

void* bolitas_new(void) {
	t_bolitas* self = (t_bolitas*) object_alloc(this_class);
	/*outlet_new(self, "signal");
	outlet_new(self, "signal");
	outlet_new(self, "signal");
	outlet_new(self, "signal");
	outlet_new(self, "signal");
	outlet_new(self, "signal");
	outlet_new(self, "signal");*/

    xo3 = outlet_new(self,"float");
    xo2 = outlet_new(self,"float");
    xo1 = outlet_new(self,"float");
    yo = outlet_new(self,"float");
    to3 = outlet_new(self,"float");
    to2 = outlet_new(self,"float");
    to1 = outlet_new(self,"float");

	random_generator.Init(1);
	random_stream.Init(&random_generator);
	note_filter.Init();
	onSampleRateChange();
	onReset();

	/*dsp_setup((t_pxobject*)self, 0);*/
	return (void *)self;
}


void bolitas_free(t_bolitas* self) {
	/*dsp_free((t_pxobject*)self);*/
}

void bolitas_dsp64(t_bolitas* self, t_object* dsp64, short* count, double samplerate, long maxvectorsize, long flags) {
	object_method_direct(void, (t_object*, t_object*, t_perfroutine64, long, void*),
						 dsp64, gensym("dsp_add64"), (t_object*)self, (t_perfroutine64)bolitas_perform64, 0, NULL);
}


void bolitas_assist(t_bolitas* self, void* unused, t_assist_function io, long index, char* string_dest) {
	if (io == ASSIST_OUTLET) {
		switch (index) {
			case 0: 
				strncpy(string_dest,"T1 Output", ASSIST_STRING_MAXSIZE); 
				break;
			case 1: 
				strncpy(string_dest,"T2 Output", ASSIST_STRING_MAXSIZE); 
				break;
			case 2: 
				strncpy(string_dest,"T3 Output", ASSIST_STRING_MAXSIZE); 
				break;
			case 3: 
				strncpy(string_dest,"Y Output", ASSIST_STRING_MAXSIZE); 
				break;
			case 4: 
				strncpy(string_dest,"X1 Output", ASSIST_STRING_MAXSIZE); 
				break;
			case 5: 
				strncpy(string_dest,"X2 Output", ASSIST_STRING_MAXSIZE); 
				break;
			case 6: 
				strncpy(string_dest,"X3 Output", ASSIST_STRING_MAXSIZE); 
				break;
		}
	}
}

void bolitas_trate(t_bolitas *x, double f)
{
	f_t_rate = f;
}

void bolitas_tclockexternalinput(t_bolitas *x, double f)
{
	f_t_clock_input_patched = f>0.5;
}

void bolitas_tclockexternal(t_bolitas *x, double f)
{
	f_t_clock_input = f>0.5;
}

void bolitas_xclockexternalinput(t_bolitas *x, double f)
{
	f_x_clock_input_patched = f>0.5;
}

void bolitas_xclocksourceinternal(t_bolitas *x, double f)
{
	x_clock_source_internal = (int) f;
}


void bolitas_xclockexternal(t_bolitas *x, double f)
{
	f_x_clock_input = f>0.5;
}

void bolitas_external(t_bolitas *x, double f)
{
	external = f > 0.5f;
}


void bolitas_tdejavu(t_bolitas *x, double f)
{
	t_deja_vu = f > 0.5f;
}

void bolitas_xdejavu(t_bolitas *x, double f)
{
	x_deja_vu = f > 0.5f;
}

void bolitas_fdejavu(t_bolitas *x, double f)
{
	f_deja_vu = constrain(f, 0.f, 1.f);
}

void bolitas_dejavulength(t_bolitas *x, double f)
{
	f_deja_vu_length = constrain(f, 0.f, 1.f);
	float deja_vu_length_index = f_deja_vu_length * (sizeof(loop_length) - 1);
	deja_vu_length = loop_length[(int) roundf(deja_vu_length_index)];
}

void bolitas_tmode(t_bolitas *x, double f)
{
	t_mode = (int) f;
}

void bolitas_xmode(t_bolitas *x, double f)
{
	x_mode = (int) f;
}

void bolitas_trange(t_bolitas *x, double f)
{
	t_range = (int) f;
}

void bolitas_xrange(t_bolitas *x, double f)
{
	x_range = (int) f;
}

void bolitas_pws(t_bolitas *x, double f)
{
	f_pulse_width_std = (float) f;
}

void bolitas_pwm(t_bolitas *x, double f)
{
	f_pulse_width_mean = (float) f;
}

void bolitas_tbias(t_bolitas *x, double f)
{
	f_t_bias = f;
}

void bolitas_xbias(t_bolitas *x, double f)
{
	f_x_bias = f;
}

void bolitas_spread(t_bolitas *x, double f)
{
	f_x_spread = f;
}

void bolitas_jitter(t_bolitas *x, double f)
{
	f_t_jitter = f;
}

void bolitas_steps(t_bolitas *x, double f)
{
	f_x_steps = f;
}


void ext_main(void* r) {
	this_class = class_new("bolitas~", (method)bolitas_new, (method)bolitas_free, sizeof(t_bolitas), NULL, A_GIMME, 0);

	class_addmethod(this_class,(method) bolitas_assist, "assist",	A_CANT,		0);
	/*class_addmethod(this_class,(method) bolitas_dsp64, "dsp64",	A_CANT,		0);*/

	class_addmethod(this_class,(method) bolitas_xrange, "xrange", A_DEFFLOAT, 0);
	class_addmethod(this_class,(method) bolitas_trange, "trange", A_DEFFLOAT, 0);
	class_addmethod(this_class,(method) bolitas_xmode, "xmode", A_DEFFLOAT, 0);
	class_addmethod(this_class,(method) bolitas_tmode, "tmode", A_DEFFLOAT, 0);
	class_addmethod(this_class,(method) bolitas_xdejavu, "xdejavu", A_DEFFLOAT, 0);
	class_addmethod(this_class,(method) bolitas_tdejavu, "tdejavu", A_DEFFLOAT, 0);
	class_addmethod(this_class,(method) bolitas_fdejavu, "fdejavu", A_DEFFLOAT, 0);
	class_addmethod(this_class,(method) bolitas_external, "external", A_DEFFLOAT, 0);
	class_addmethod(this_class,(method) bolitas_xclockexternal, "xclockexternal", A_DEFFLOAT, 0);
	class_addmethod(this_class,(method) bolitas_tclockexternal, "tclockexternal", A_DEFFLOAT, 0);
	class_addmethod(this_class,(method) bolitas_tclockexternalinput, "tclockexternalinput", A_DEFFLOAT, 0);
	class_addmethod(this_class,(method) bolitas_xclockexternalinput, "xclockexternalinput", A_DEFFLOAT, 0);
	class_addmethod(this_class,(method) bolitas_trate, "trate", A_DEFFLOAT, 0);
	class_addmethod(this_class,(method) bolitas_xbias, "xbias", A_DEFFLOAT, 0);
	class_addmethod(this_class,(method) bolitas_tbias, "tbias", A_DEFFLOAT, 0);
	class_addmethod(this_class,(method) bolitas_spread, "spread", A_DEFFLOAT, 0);
	class_addmethod(this_class,(method) bolitas_steps, "steps", A_DEFFLOAT, 0);
	class_addmethod(this_class,(method) bolitas_jitter, "jitter", A_DEFFLOAT, 0);
	class_addmethod(this_class,(method) bolitas_xclocksourceinternal, "xclocksourceinternal", A_DEFFLOAT, 0);
	class_addmethod(this_class,(method) bolitas_dejavulength, "dejavulength", A_DEFFLOAT, 0);
	class_addmethod(this_class,(method) bolitas_pwm, "pwm", A_DEFFLOAT, 0);
	class_addmethod(this_class,(method) bolitas_pws, "pws", A_DEFFLOAT, 0);
	class_addmethod(this_class,(method) bolitas_run, "run", A_DEFFLOAT, 0);

	/*class_dspinit(this_class);*/
	class_register(CLASS_BOX, this_class);
}



