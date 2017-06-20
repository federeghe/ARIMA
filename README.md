# ARIMA Model Estimation

This software is derived from TISEAN suite and it enables the possibility to estimate a multi-variate
ARIMA model and forecast at 1-step the following values.

### Configuration

Check the file `config.h`.
Most of the options are configurable only at compile time to performance reason. Pay attention to the
memory consumption if you are in a low-memory system (e.g. embedded system).

The verbosity for debugging purposes can be changed with the `ARIMA_DEBUG` and `ARIMA_DEBUG_EXTRA` options.



### Build&Run instructions

In order to build the example application configure the Makefile with your compiler and linker.
Then simply run:

`make`

And then you can run the application (you have to provide `temp1.txt` and `temp2.txt` files)


### License

GNU General Public License 2.0

### Contributors

* Authors of [TISEAN suite](https://www.pks.mpg.de/~tisean/)
* [Federico Reghenzani](https://home.deib.polimi.it/reghenzani)
