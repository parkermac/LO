# LiveOcean Forecast Operators Manual

#### This is an introduction to the LiveOcean daily forecast, how to check on it every morning, and what to do when something goes wrong. It is written assuming you already know about our analysis machines (apogee and perigee) and hyak (mox and klone).

The LiveOcean forecast has been running continuously since late 2015. Its results are used by a half-dozen stakeholders who rely on it for anything from public shellfish safety decisions to open boundary conditions for other forecast models. In all this time we have only fully missed a few forecast days, and forecasts have only been significantly late a few dozen times. Maintaining a reliable forecast is essential. We continuously improve the reliability of the system, and we check on the forecast every morning.

---

## Morning Ritual

If you are on the list of operators you will receive an email every morning shortly after 7 AM local time titled: **[klone-hyak] LO forecast klone**. The contents will look something like this:
```
Running ROMS forecast 2023.08.03-2023.08.05

======== f2023.08.03 =========
 > started at 2023.08.03 03:00:02
/gscratch/macc/parker/LO_roms/cas6_traps2_x2b/f2023.08.03
 - time to get forcing = 259 sec
========= Finished? ==========
----------- stderr -----------
slurm_load_jobs error: Unable to contact slurm controller (connect failure)

 - time to run ROMS = 2896 sec
 - ROMS SUCCESS
 - copy & clean up = 197 sec
 > finished at 2023.08.03 03:56:00

======== f2023.08.04 =========
 > started at 2023.08.03 03:56:00
/gscratch/macc/parker/LO_roms/cas6_traps2_x2b/f2023.08.04
 - time to get forcing = 41 sec
 - time to run ROMS = 841 sec
 - Run blew up, blow ups = 0
 - blew up at 2023.08.03 04:10:48
 - time to run ROMS = 3247 sec
 - ROMS SUCCESS
 - copy & clean up = 195 sec
 > finished at 2023.08.03 05:08:17

======== f2023.08.05 =========
 > started at 2023.08.03 05:08:17
/gscratch/macc/parker/LO_roms/cas6_traps2_x2b/f2023.08.05
 - time to get forcing = 44 sec
 - time to run ROMS = 2896 sec
 - ROMS SUCCESS
 - copy & clean up = 209 sec
 > finished at 2023.08.03 06:00:53
```
If you see **SUCCESS** on each of the three days, **then you are done**. The run worked fine, and it is unlikely that there will be any problems with the post-processing. You can see a couple of hiccups in this log file. On the first day it reported a slurm_load_jobs error. This happens occasionally, and the python driver on klone just keeps trying until it works. You can also see that on the second day that the run blew up. This is common during strong spring tides, as was the case here, and the python driver on klone just runs the day again with a smaller time step.

#### If you do not see SUCCESS three times in the email, then there was some problem. Here are the possibilities and what to do about them.

(1) **It may be a klone problem.** Occasionally it will just not start the job and then eventually the python driver gives up. However, we try to run the forecast THREE TIMES every day, at 3 AM, 8 AM, and 5 PM. It is highly like that it will work on the second try.

(2) **It may be a forcing problem.** If this is the case you will see a lot of error messages in the log file about forcing files. If this happens, logon to apogee and go to /dat1/parker/LO/driver, and do  `cat *new.log` which will produce something like:
```
  frc=atm00, day=2023.08.03, result=success, note=NONE
  start=2023.08.03 02:55:01 (took 89 sec)
  /dat1/parker/LO_output/forcing/cas6/f2023.08.03/atm00

  frc=ocn00, day=2023.08.03, result=success, note=planA
  start=2023.08.03 01:20:03 (took 523 sec)
  /dat1/parker/LO_output/forcing/cas6/f2023.08.03/ocn00

  frc=tide00, day=2023.08.03, result=success, note=NONE
  start=2023.08.03 02:00:02 (took 16 sec)
  /dat1/parker/LO_output/forcing/cas6/f2023.08.03/tide00

  frc=traps2, day=2023.08.03, result=success, note=NONE
  start=2023.08.03 02:05:02 (took 307 sec)
  /dat1/parker/LO_output/forcing/cas6/f2023.08.03/traps2
```
These are the log files generated every morning on apogee by `LO/driver/driver_forcing3.py`. In this example everything worked (you see "result=success" for all four). If one fails ("result=fail") it may give you a clue in the "note". The default note is "NONE" except for ocn00 where the default is "planA". **ocn00** can also generate "planB" results where it just used a different method to download the HYCOM files, or "planC" where is copies over the previous day's ocn00 forcing and adds a day to the last time index. All three will work to run the model, but if planC happens more than a few days in a row then to ocean boundary condition is not current. **atm00** can also generate "planB" which uses the forcing from the previous day and persists the third day. If there is an atm00 planB you should notify David Darr and he can check on our source for the WRF files.

Because the system is designed to be reliable, in most cases even if it used planB or planC you will still see "result=success" and the forcing is fine to run the model. If there was a "result=fail" you will have to re-run the forcing by hand on apogee. Because of permission issues I think only David Darr can do this. An example of a command to run, from /dat1/parker/LO/driver on apogee is:
```
python driver_forcing3.py -g cas6 -r forecast -f ocn00 > ocn00_new.log &
```
_**David: we need to check that you can do this.**_

Once the missing forcing is generated you can just wait for klone to try the forecast again.

(3) **A final thing to check on is the post-processing**, but this is optional. You would want to do it if you noticed that the website movies were not updating (e.g. look at the datestamp in the lower left corner of the movie http://faculty.washington.edu/pmacc/LO/p5_full_salt_top.html). Logon to apogee and go to /dat1/parker/LO/driver, and do  `cat post1.log.log` which will produce something like:
```
---------- Post-processing forecast for 2023.08.03----------
All files found. Beginning post-processing.

nest_wgh 2023.08.03 SUCCESS
start=2023.08.03 06:02:06 (took 542 sec)
/dat1/parker/LO_output/post/cas6_traps2_x2b/f2023.08.03/nest_wgh

surface1 2023.08.03 SUCCESS
start=2023.08.03 06:11:09 (took 166 sec)
/dat1/parker/LO_output/post/cas6_traps2_x2b/f2023.08.03/surface1

layers1 2023.08.03 SUCCESS
start=2023.08.03 06:13:55 (took 610 sec)
/dat1/parker/LO_output/post/cas6_traps2_x2b/f2023.08.03/layers1

ubc1 2023.08.03 SUCCESS
start=2023.08.03 06:24:06 (took 168 sec)
/dat1/parker/LO_output/post/cas6_traps2_x2b/f2023.08.03/ubc1

sequim1 2023.08.03 SUCCESS
start=2023.08.03 06:26:56 (took 112 sec)
/dat1/parker/LO_output/post/cas6_traps2_x2b/f2023.08.03/sequim1

critfc1 2023.08.03 SUCCESS
start=2023.08.03 06:28:48 (took 1174 sec)
/dat1/parker/LO_output/post/cas6_traps2_x2b/f2023.08.03/critfc1

daymovie0 2023.08.03 SUCCESS
start=2023.08.03 06:48:23 (took 356 sec)
/dat1/parker/LO_output/post/cas6_traps2_x2b/f2023.08.03/daymovie0

drifters0 2023.08.03 SUCCESS
start=2023.08.03 06:54:20 (took 893 sec)
/dat1/parker/LO_output/post/cas6_traps2_x2b/f2023.08.03/drifters0

archive0 2023.08.03 SUCCESS
start=2023.08.03 07:09:14 (took 716 sec)
/dat1/parker/LO_output/post/cas6_traps2_x2b/f2023.08.03/archive0

Total time for all post jobs = 4744.9 sec
```
The code that generated this log file is `LO/driver/driver_post1.py`. It will keep looking for the forecast files for 15 hours (until 8 PM if started at 5 AM), and only start when the last of the history files shows up.

(4) **If all else fails**, contact David Darr at darrd@uw.edu.

(5) **Don't despair!** Even if the forecast completely fails _two days in a row_, it will still be able to run on the third day using the last good restart file.

---

#### Background and Extra Info

**Crontabs:** Look in /LO/driver/crontabs for the crontabs currently running on apogee and klone.

Here current daily lineup of the forecasts and post-processing. Times in parentheses () only run if the previous forecast failed.

| Job | time | time | time | time | time | time |
| --- | --- | --- | --- | --- | --- | --- |
| cas6_traps2_x2b | 3 AM | | (8 AM) | | | (5 PM) |
| post1 | | 5 AM | | | | |
| wgh2_t0_xnb | | | | 1 PM | | |
| post2 | | | | | 2 PM | |


**How does the driver know if the forecast already ran on klone?** Since we re-run the python driver three times every day on klone using cron, there is a little file generated by the driver when the forecast completes successfully, e.g. `LO/driver/forecast3_done_2023.08.03.txt`. If the driver finds this then it exits.

**People who use the forecast products:** In the event that the forecast will be later than about Noon local time, send a brief email with an expected forecast time to this list (just paste this text string into the "To" section of the email):
```
adioso@uw.edu;evan@scootscience.com;djlatornell@gmail.com;Zhaoqing.Yang@pnnl.gov;ryan.mccabe@noaa.gov;harpers@uw.edu;pmb@uw.edu;p.maccready@gmail.com
```
You can see more info about the post processing and the people who use it in LO/notes/post_users.md.

**Hyak maintenance:** the hyak system has scheduled maintenance on the second Tuesday of every month. They will stop or kill jobs that have scheduled run times that go beyond 9 AM. We have worked hard to make sure that, even when there are multiple blow-ups, the last forecast day will still start before 7 AM. Then because `#SBATCH --time=02:00:00` in its `LO/driver/batch/klone3_batch_BLANK.sh` it will finish before they shut it is down. One of the more aggravating issues is when a forecast fails on a maintenance day. Then you have to wait until the late afternoon before you can rerun ROMS by hand, or wait for the last crontab attempt.

---
