# Energy Weighted Observable Correlations (EWOCs)

[![standard-readme compliant](https://img.shields.io/badge/standard--readme-OK-green.svg?style=flat-square)](https://github.com/RichardLitt/standard-readme)

A set of C++ and python-based tools for calculating energy-weighted differential cross sections involving pairwise jet and subjet observables; I playfully call these differential cross sections/correlations "Energy Weighted Observable Correlations" (EWOCs).


## Table of Contents

- [Usage](#usage)
  - [Setup](#setup)
  - [Visualization](#vis)

- [File Structure](#structure)

- [Physics Background](#background)
  - [Jets](#jets)
  - [Sub-jets](#subjets)

- [Maintainers](#maintainers)
- [Contributing](#contributing)
- [License](#license)


## <a name="usage"></a> Usage

### <a name="setup"></a> Setup

To set up, go to the ``Makefile.inc`` file and edit the locations of the dependencies to ensure they are accurate for your system. Then, it should be enough to use the  ``make`` command of the main directory of this project. This should prepare the executables and scripts that I have been using to understand EWOCS.

I note that I follow the unconventional choice of using the fish shell language in most of these scripts; sorry for any inconvenience.


### <a name="vis"></a> Visualization
To visualize an example event, try ``./scripts/vis_event_single.sh --process qcd --level parton`` from the main directory. Other options can also be provided; look at this script and the functions it calls to find out more: ``cat scripts/vis_events_single.sh``.


## <a name="structure"></a> File Structure

I think it might be helpful even for myself to have some notes regarding the file structure here; I hope if you're visiting that it helps you too!

The files in this repository are usually grouped as useful either for "writing", in C++, (data analysis and storage) or "plotting", in python (data visualization). These are stored in the ``write_tools`` folder and the ``plot_tools`` folder, respectively:

<p align="center">
<img src="./doc_figs/basic_structure.png" width="500">
</p>

The writing and plotting utilities are grouped together into the executables ``write_tools/write_ewocs`` and ``plot_tools/plot_ewocs``, respectively. Events can also be visualized with the executable ``plot_tools/event_vis``. These can all be prepared using the ``make`` command in the main folder:

<p align="center">
<img src="./doc_figs/make.png" width="300">
</p>

Finally, the scripts in the ``scripts`` folder tie the writing and plotting tools, so that one can more easily understand EWOCs, with fewer presses.

<p align="center">
<img src="./doc_figs/scripts.png" width="800">
</p>

I note (again) that I follow the unconventional choice of using the fish shell language in most of these scripts; sorry for any inconvenience.


## <a name="background"></a> Physics Background
Energy weighting is a way of taking the physics at certain energy scales "more seriously": observables which include a linear energy weighting give higher energy particles more weight than they do to low energy particles.
For this reason, energy weighting is useful in situations when we don't understand dynamics at particular energy scales, and want to isolate physics at energy scales we _do_ understand.

In the case of Quantum Chromodynamics (QCD), we want to test the behavior of the "fundamental" building blocks of QCD at high energy scales: quarks and gluons. Energy weighting allows us to isolate the high energy physics of quarks and gluons while remaining more agnostic to the behavior of QCD at low energy scales -- a difficult problem which has evaded physicists for a long time, and which we would like to avoid for cleanness and simplicity.

Energy weighted observables are useful in QCD, and have been for a long while. In part, this is because they are more easily calculable than other observables in QCD; in part again, this is because they are less susceptible to complicated contributions from soft radiation: contributions of soft radiation to observables is naturally suppressed by a linear energy weighting (or an energy weighting of any positive power)!

### <a name="jets"></a> Jets
Quarks and gluons are difficult to understand both theoretically and experimentally; the concept of the **jet** was introduced to enable further understanding of these complicated objects.

Due to a phenomenon called "asymptotic freedom", quarks and gluons are easier to understand at high energies; they behave almost as free, non-interacting particles, streaming by one another with nary a care. Of course, we do not see them in our daily lives; this is due to a phenomenon called "confinement" -- a rough statement that we still don't fully understand after 70 years, and which states that quarks and gluons interact more and more strongly at lower and lower energies, until they finally get bound up into other, composite particles such as protons and neutrons.

Jets help us bridge the gap between the low energy reality of our world and the high energy fantasy of free quarks and gluons. Roughly, the idea is as follows:
    1) To understand the high energy structure of our universe, we slam particles together at high energies and try to understand what comes out.
    2) Quarks and gluons "appear" (read: "are effective descriptions of what we see") at high energies, and our theory of quarks and gluons may be used to describe physical phenomena for a very, VERY, short time after we smash our collidees together.
    3) As time passes, quarks and gluons become less effective descriptions of our system. Again, for a short time, the quarks and gluons "spread out", and they radiate away energy through other quarks and gluons. I believe we do not yet have a full description of this phenomena, but it is _somewhat_ similar to the way an quickly charged particle radiates energy in the form of an external electric field. A quickly moving electron releases photons; analogously, a quickly moving quark or gluon releases gluons and other forms of QCD radiation.
    4) The structures that emerge in this radiative process are roughly cone-shaped -- because the "initial quark or gluon" is moving very fast, the radiation it releases tends to be moving very fast in the same direction.
    5) A **jet algorithm** is, loosely, an algorithm which groups together particles which lie together in cones. Jet algorithms are usually designed so that they can effectively represent the "original quarks and gluons" of step 2. When we apply a jet algorithm on a set of particles coming out of a particular particle collision, we call the resulting group of particles a **jet**. We can use our theory of quarks and gluons, and some complicated physics technology, to accurately describe the properties of jets that we see in particle collisions.

Jet algorithms are usually defined with a _radius_ -- roughly, they search for particles which lie in a circle in the space of outgoing particle angles.

Here is an visualization of jets of radius 1.5 in an example particle collision, generated by calling ``./scripts/vis_event_scattter.sh -l hadron -p qcd -j akt -R 1.5 -r 0.0``:
![plot](./doc_figs/example_jet_vis.png)

In this picture, higher energy jets are more red, and lower energy jets are more blue. The points indicate particles in the jet and, roughly, the light background indicates in which jet the jet algorithm would place an additional, low energy particle in the event.


Examining energy-weighted observables within jets is a popular and intriguing way to try and understand the energy profile of the "initial quark or gluon", by understanding what we can about correlations in this radiation profile.


### <a name="jets"></a> Sub-jets

Sub-jets are jets within jets. When we consider all the particles within a jet of radius R (i.e. grouped together with a jet algorithm of radius R) and then apply _another_ jet algorithm with radius r < R, we are left with **subjets**. Energy weighted correlations of pairwise observables of subjets within a jet (subjet energy weighted observable correlations, or subjet EWOCs) provide an unexplored way of probing the patterns of energy that emerge in QCD.

Here is an example visualization of some sub-jets of the highest energy jet from above (the dark red jet), generated by calling ``./scripts/vis_event_scatter.sh -l hadron -p qcd -j akt -R 2 -r 0.5``:
![plot](./doc_figs/example_subjet_vis.png)

## <a name="maintainers"></a> Maintainers

[@samcaf](https://github.com/samcaf)


## <a name="contributing"></a> Contributing

PRs accepted.

Small note: If editing the README, please conform to the [standard-readme](https://github.com/RichardLitt/standard-readme) specification.


## <a name="license"></a> License

MIT © 2021 Samuel Alipour-fard
