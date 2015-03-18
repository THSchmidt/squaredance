![http://squaredance.googlecode.com/files/squaredance_logo3.png](http://squaredance.googlecode.com/files/squaredance_logo3.png)

# Cross-sectional Area Calculation #

**squaredance** is a grid-based approach to calculate 2D-slice-wise cross-sectional areas of proteins or any other group of atoms defined via GROMACS NDX groups.

## Features ##
  * Grid-based calculation of...
    * ...molecular cross-sectional areas, considering the atoms' van der Waals radii, per slice
    * ...the total volume of internal cavities per slice
    * ...the molecular surface (vdW, SAS) per slice
  * Projecting atoms within a defined range to a 2D plane, e.g. for calculating the area occupied by the selected atoms


## News ##

_2012-09-xx_ The first release of squaredance is available for download.


## Quick start ##
  * Download and unpack the squaredance archive
  * Copy the folder to your lib directory (e.g. `cp -r squaredance/ ~/lib`)
  * Make squaredance executable (e.g. `chmod 755 ~/lib/squaredance/squaredance`)
  * Set a soft link to the squaredance executable (e.g. `ln -s ~/lib/squaredance/squaredance ~/bin/squaredance`)

You can run squaredance by typing
```bash
squaredance -f protein.gro```

Get the program help by typing
```bash
squaredance -h```


## Documentation ##

We provide **[installation instructions](Installation.md)** and a **[manual](Manual.md)** on this program's project page.


## References ##
If you use squaredance for **publication**, please cite

> http://code.google.com/p/squaredance


**Electronic documents** should include a direct link to the squaredance hosting page:

> http://code.google.com/p/squaredance


For **posters** or nerds we provide a QR Code, referencing to this Google Project Page:

  * [QR Code 96dpi](http://code.google.com/p/squaredance/downloads/detail?name=squaredance_qr96dpi.png)
  * [QR Code 300dpi](http://code.google.com/p/squaredance/downloads/detail?name=squaredance_qr300dpi.png)


## Related projects ##

  * http://code.google.com/p/lambada-align
  * http://code.google.com/p/inflategro2