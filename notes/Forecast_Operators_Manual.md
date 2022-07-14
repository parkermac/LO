# LiveOcean Forecast Operators Manual

#### This is an introduction to the LiveOcean daily forecast, how to check on it every morning, and what to do when something goes wrong. It is written assuming you already know about our analysis machines (apogee and perigee) and hyak (mox and klone).

---

## Morning Ritual


The LiveOcean forecast has been running continuously since late 2015. Its results are used by a half-dozen stakeholders who rely on it for anything from public shellfish safety decisions to open boundary conditions for other forecast models. In all this time we have only fully missed a few forecast days, and forecasts have only been significantly late a few dozen times. Maintaining a reliable forecast is essential. We continuously improve the reliability of the system, and we check on the forecast every morning.

_**Here are the typical daily steps:**_

- First, open a browser and look at https://liveocean.apl.uw.edu/output/. Anytime after about 6 AM you should see a folder at the bottom of the list named after the current day f[YYYY.MM.DD] (e.g. f2022.07.12). This is the public server where we make daily forecast post-processing files available to users. **If the folder is there, everything probably worked fine and if you are in a hurry you can stop here.**
- At 6:30 AM official operators should get two emails. One, titled "LOcean results", is the screen output of driver_roms1.py from the primary forecast run by Parker on klone. The other, titled "LO_mox", is from the backup forecast run by Parker on mox. These give information on each of the three forecast days. For example, the lower third of "LOcean results" today has info about day three of the forecast:
```
OOOOOOOOOOOOOOOOOOOOOOO f2022.07.14 OOOOOOOOOOOOOOOOOOOOOOOO
 > started at 2022.07.12 04:52:14
 - force_dir:    /gscratch/macc/parker/LO_output/forcing/cas6_v0/f2022.07.14
 - dot_in_dir:   /gscratch/macc/parker/LO/dot_in/cas6_v0_u0kb
 - dot_in_shared_dir: /gscratch/macc/parker/LO/dot_in/shared
 - roms_out_dir: /gscratch/macc/parker/LO_roms/cas6_v0_u0kb/f2022.07.14
 - log_file:     /gscratch/macc/parker/LO_roms/cas6_v0_u0kb/f2022.07.14/log.txt
 - roms_ex_dir:  /gscratch/macc/parker/LiveOcean_roms/makefiles/u0kb
==================== Copy forcing atm0 =====================
==================== Copy forcing ocn0 =====================
==================== Copy forcing riv0 =====================
==================== Copy forcing tide0 ====================
 - time to get forcing = 58 sec
 - Blow-ups = 0
======================= Make dot in ========================
=================== Create batch script ====================
========================= Run ROMS =========================
 - time to run ROMS = 170 sec
 - llcount = 0
 - log done and closed
 - Run blew up, blow ups = 0
 - blew up at 2022.07.12 04:56:07
 - log.txt file saved to /gscratch/macc/parker/LO_roms/cas6_v0_u0kb/f2022.07.14_blowup
 - ocean_his_0001.nc saved to /gscratch/macc/parker/LO_roms/cas6_v0_u0kb/f2022.07.14_blowup
 - Blow-ups = 1
======================= Make dot in ========================
=================== Create batch script ====================
========================= Run ROMS =========================
 - time to run ROMS = 2541 sec
 - llcount = 0
 - log done and closed
 - ROMS SUCCESS
=== Make output directory on apogee.ocean.washington.edu ===
===== Copy ROMS output to apogee.ocean.washington.edu ======
 - time to move history files and clean up = 187 sec
 > finished at 2022.07.12 05:41:39
```
This is from a strong spring tide day, and the model often blows up once or even twice during springs. Things I watch for are excessive "time to get forcing" (e.g. over 300 sec) which may indicate that klone is having some problem. If there is a problem, send an email describing it to help@uw.edu with "hyak" in the subject line.
- To be thorough, I also check on the forcing to see if any "planB" actions were taken. To do this, logon to apogee and go to `/dat1/parker/LO/driver` and do `cat *0.log` which will produce something like:
```
  frc=atm0, day=2022.07.12, result=success, note=NONE
  start=2022.07.12 02:55:02 (took 61 sec)
  /dat1/parker/LO_output/forcing/cas6_v0/f2022.07.12/atm0

  frc=ocn0, day=2022.07.12, result=success, note=planA
  start=2022.07.12 01:20:07 (took 590 sec)
  /dat1/parker/LO_output/forcing/cas6_v0/f2022.07.12/ocn0

  frc=riv0, day=2022.07.12, result=success, note=NONE
  start=2022.07.12 02:05:02 (took 355 sec)
  /dat1/parker/LO_output/forcing/cas6_v0/f2022.07.12/riv0

  frc=tide0, day=2022.07.12, result=success, note=NONE
  start=2022.07.12 02:00:02 (took 24 sec)
  /dat1/parker/LO_output/forcing/cas6_v0/f2022.07.12/tide0
```
These are from the four times `driver_forcing.py` was run in the early morning by cron. In this case everything looks good. The default for ocn0 is that the note says planA, and for all others the default is NONE. Both ocn0 and atm0 have a planB, and ocn0 even has a planC. If planB is activated for atm0 that means there was a problem getting the WRF files and this means you should alert David Darr at darrd@uw.edu. If planC is activated for ocn0 that means it had a problem getting HYCOM files, and you should rerun it by hand so that we have a better fallback for the next day, even if it has "result=success".
- If any of the forcings failed ("result=fail" or just missing log results) then the forecast also almost certainly did not run. In this case the standard procedure is to:
 - Look at /dat1/parker/LO_output/forcing/cas6_v0/f[YYYY.MM.DD]/[frc]0/Info/screen_output.txt for a clue about what the problem may have been.
 - Often you can just rerun the forcing on apogee by hand with a command like:
 - `python driver_forcing.py -g cas6 -t v0 -r forecast -s continuation -f atm0`
 - Then logon to klone and first check that nothing else is running that would gum things up (use top and squeue). Then rerun the forecast by hand with a command like:
 - `python3 driver_roms1.py -g cas6 -t v0 -x u0k -r forecast -np 400 -N 40 < /dev/null > rerun_[YYYY.MM.DD].log &`
 - Assuming the forecast finishes by 1 PM the automated `post_process1.py` on apogee should work on its own. If it is later than 1 PM you will have to run this by hand as well with a command like:
 - `python driver_post1.py -gtx cas6_v0_u0kb -r forecast -ro 0 < /dev/null > post_rerun_[YYYY.MM.DD].log &`
 - If the primary forecast failed but the **backup** worked you can first kill any running instance of driver_post1 on apogee, and then run the post_processing on the backup using a command like:
 - `python driver_post1.py -gtx cas6_v0_u0mb -r forecast -ro 0 < /dev/null > post_rerun_[YYYY.MM.DD].log &`

_**These steps will deal with about 95% of all error modes!**_

---

## Background info on how the forecast system works

Every morning a series of cron jobs run on apogee, klone, and mox in Parker's account that make the LiveOcean three-day forecast.  This spans today, tomorrow, and the day after, and like all LO ROMS runs, it is done in three separate one-day jobs on hyak. The system consists of three categories:
- Preprocessing
- ROMS
- Post-processing

---

#### Preprocessing

This is creating the forcing for the run.  The current lineup is tide0, riv0, ocn0, and atm0, and these are run in Parker's account on apogee by the `driver_forcing.py` lines in his crontab:

```
LOd=/dat1/parker/LO/driver
HOSTNAME=apogee
# * * * * * source ~/.bashrc; python $LOd/test_Ldir.py > ~/test.log
# * * * * * date > ~/date.log
00 02 * * * source ~/.bashrc; python $LOd/driver_forcing.py -g cas6 -t v0 -r forecast -f tide0 > $LOd/tide0.log
05 02 * * * source ~/.bashrc; python $LOd/driver_forcing.py -g cas6 -t v0 -r forecast -f riv0 > $LOd/riv0.log
20 01 * * * source ~/.bashrc; python $LOd/driver_forcing.py -g cas6 -t v0 -r forecast -f ocn0 > $LOd/ocn0.log
55 02 * * * source ~/.bashrc; python $LOd/driver_forcing.py -g cas6 -t v0 -r forecast -f atm0 > $LOd/atm0.log
00 05 * * * source ~/.bashrc; python $LOd/driver_post1.py -gtx cas6_v0_u0kb -ro 0 -r forecast < /dev/null > $LOd/post1.log
```

The numbers at the start of each line are MM HH (minute and hour), so you can see that these are all run between 1:20 AM and 2:55 AM local time (local time). We run atm0 last because the job (handled by David Darr) to get WRF files from UW Atmospheric Sciences is the last to finish.

---

#### ROMS

_Under construction..._

One issue we need to deal with is what non-Parker users have to do to run the forecast.

NOTE: the hyak system has scheduled maintenance on the second Tuesday of every month. They will stop or kill jobs that have scheduled run times that go beyond 9 AM. We have worked hard to make sure that, even when there are multiple blow-ups, the last forecast day will still start before 7 AM. Then because `#SBATCH --time=02:00:00` in its `LO/driver/batch/[klone,mox]1_batch_BLANK.sh` it will finish before they shut it is down. One of the more aggravating issues is when a forecast fails on a maintenance day. Then you have to wait until the late afternoon before you can rerun ROMS by hand.

---

#### Post-processing

This is orchestrated by the last line in the crontab above by `driver_post1.py`. This makes all the files that end up on the LiveOcean server at APL https://liveocean.apl.uw.edu/output/. It also makes movies that are pushed to Parker's LiveOcean website http://faculty.washington.edu/pmacc/LO/LiveOcean.html.

Issues to deal with will be:
- Do other users have permission to scp files to the APL server?
- Do other users have permission to scp files to homer (where Parker's website files are)?

If you have to rerun ROMS late in the day, you can run `driver_post1.py` on apogee at the same time - it will keep looking for 8 hours, and only start when the last of the history files shows up.