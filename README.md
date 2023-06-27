# Guitar synthesis using physical modelling

This repository contains two scripts:
1. `ks_acoustic_guitar.m` — acoustic guitar synthesis using Karplus-Strong algorithm.
3. `fdtd_acoustic_guitar.m` — acoustic guitar synthesis using finite difference method (including string collisions with fretboard and finger).

Since simulation times for finite difference method can be long, two demo recordings are included:
1. `fdtd_acoustic_guitar_open_strings.wav` — six open strings are exited by a pluck of increasing amplitude (in the end you can hear string collisions with fretboard).
2. `fdtd_acoustic_guitar_chord_tapping.wav` — six string are tapped against the fretboard by a finger to produce a Fmaj7 chord (you can hear rattling from finger/string/fretboard collisions).

All of the relevant references are provided within MATLAB scripts.
