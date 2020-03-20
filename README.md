* **Input** an original probe image

* **output** SH coefficients **&** the reconstructed image using these SH coefficients.

*You can also change `SH_ORDER` in SH.h*



```bash
$ ./SHTools.exe [options] <input>
<input>              .hdr probe image
Options: --help      Print help message
         -h          Hanning  filter
         -l          Lanczos  filter
         -g          Gaussian filter
```
