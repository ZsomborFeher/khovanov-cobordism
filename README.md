# Khovanov Cobordism Calculator

A Python module to calculate cobordism maps induced on Khovanov homology.

Written by Zsombor Fehér, 2024.

## Quick Start Guide

Example: Distinguishing two ribbon disks of 6_1, given as two band diagrams on the same link diagram.

1.  Make sure you have SciPy and NumPy installed and updated to the latest version:
    ```python
    pip install --upgrade scipy
    ```
    (Some methods require [Sage](https://www.sagemath.org/) to be installed, but tasks like this example do not.)

2.  Start a Python session, copy the file `khovanov.py` to your working folder, and import that module:
    ```python
    from khovanov import *
    ```

3.  Obtain a PD code of the link, and determine how the PD code corresponds to the diagram.

    Label each crossing in your diagram with its corresponding index in the PD code.
    Mark the 0th strand of each crossing with a dot, corresponding to the 0th element of its PD code entry.
    (Since this module does not have a graphical interface yet, you might find it useful to draw the link in
    SnapPy's Plink Editor, and enable Info -> DT labels and Info -> PD code for getting this correspondence.)

4.  Create a `Link` object using the PD code:
    ```python
    L = Link([(9, 4, 10, 5), (5, 8, 6, 9), (11, 2, 12, 3), (3, 10, 4, 11), (1, 7, 2, 6), (7, 1, 8, 12)])
    ```
    You can verify that your diagram is marked correctly using `print(L)`.

5.  Create the two `Cobordism` objects based on the band diagrams and your marking of crossings.

    For example, `S0.band_move(-1, (0, 0), (2, 1))` represents a (-1)-twisted band, connecting the 0th crossing's
    0th strand to the 2nd crossing's 1st strand, positioned on the band's left side.
    ```python
    S0 = Cobordism(L)
    S0.band_move(-1, (0, 0), (2, 1))
    S0.finish()
    S1 = Cobordism(L)
    S1.band_move(-1, (1, 2), (3, 3))
    S1.finish()
    ```
    You can check the resulting movie of the cobordism with `print(S0)`.

6.  Calculate the Khovanov-Jacobsson classes (`CKhElement` objects) and compare them in homology:
    ```python
    compare(S0.KJ_class(), S1.KJ_class())
    ```

7.  If this fails to distinguish the surfaces, try mirroring the cobordisms:
    ```python
    compare(S0.mirror().KJ_class(), S1.mirror().KJ_class())
    ```

See the top of the source code for more details about classes and their methods.
