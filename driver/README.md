# README for driver

## This code is the drivers for creating forcing, running ROMS, and other jobs that require us to control complex, repetitive tasks, either for a single forecast or for a long series of days.

---

### `driver_roms_mox.py`

This runs ROMS for a single forecast or for many days. It is organized to use the LO run naming system: [gtagex] = [gridname]\_[tag]\_[ex_name].

NOTE: This is much improved from the LiveOcean version:
- Being implemented in python instead of a shell script makes the code easier to follow.
- There are a number of optional flags to allow testing of each of the elements.
- The screen output is cleaned up, and there is a "verbose" mode.
- We use Path objects.
- There is a new [tag_alt] logic which allows you to add new forcing variations to some existing [gtag], but do all the ROMS writing to a different [gtag]. The [tag_alt] flag refers to the existing one where the forcing is read from.

NOTE: This code hardwires the fact that it is being run on mox and getting its forcing from boiler in /data1/parker/LiveOcean_output. I wonder if anyone else can run the scp steps that use "parker@boiler"?

---

### `driver_forcing.py`

This runs any of the forcing jobs, for one or more days, for any [gtag].
