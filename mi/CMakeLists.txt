###############################
cmake_minimum_required(VERSION 3.12)

# MEC library
project(MIPlaitsLib)

set(PLAITS_SRC
        plaits/dsp/engine/additive_engine.cc
        plaits/dsp/engine/bass_drum_engine.cc
        plaits/dsp/engine/chord_engine.cc
        plaits/dsp/engine/fm_engine.cc
        plaits/dsp/engine/grain_engine.cc
        plaits/dsp/engine/hi_hat_engine.cc
        plaits/dsp/speech/lpc_speech_synth.cc
        plaits/dsp/speech/lpc_speech_synth_controller.cc
        plaits/dsp/speech/lpc_speech_synth_phonemes.cc
        plaits/dsp/speech/lpc_speech_synth_words.cc
        plaits/dsp/engine/modal_engine.cc
        plaits/dsp/physical_modelling/modal_voice.cc
        plaits/dsp/speech/naive_speech_synth.cc
        plaits/dsp/engine/noise_engine.cc
        plaits/dsp/engine/particle_engine.cc
        plaits/dsp/physical_modelling/resonator.cc
        plaits/resources.cc
        plaits/dsp/speech/sam_speech_synth.cc
        plaits/dsp/engine/snare_drum_engine.cc
        plaits/dsp/engine/speech_engine.cc
        plaits/dsp/string.cc
        plaits/dsp/engine/string_engine.cc
        plaits/dsp/physical_modelling/string_voice.cc
        plaits/dsp/physical_modelling/string.cc
        plaits/dsp/engine/swarm_engine.cc
        plaits/dsp/engine/virtual_analog_engine.cc
        plaits/dsp/voice.cc
        plaits/dsp/engine/waveshaping_engine.cc
        plaits/dsp/engine/wavetable_engine.cc
)

set(MI_COMMON_SRC
       stmlib/dsp/units.cc
       stmlib/dsp/atan.cc 
       stmlib/utils/random.cc
)

set(MILIB_SRC 
        ${MI_COMMON_SRC}
        ${PLAITS_SRC}
        )

include_directories(
        "${PROJECT_SOURCE_DIR}"
)

#option(BUILD_SHARED_LIBS "" ON)
add_library(MIPlaitsLib ${MILIB_SRC} )
#GENERATE_EXPORT_HEADER(
#    MUTABLELib  # generates the export header shared_EXPORTS.h automatically
#    BASE_NAME MUTABLELib
#    EXPORT_MACRO_NAME MUTABLELib_EXPORTS
#    EXPORT_FILE_NAME MUTABLELib_EXPORTS.h
#    STATIC_DEFINE MUTABLELib_BUILT_AS_STATIC)
#add_library(MUTABLELib ${MILIB_SRC} )

