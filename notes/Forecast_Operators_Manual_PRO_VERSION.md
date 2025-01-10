# README for Running the LiveOcean Forecast Model

## These notes were started 2024.12.07 by Parker MacCready, and are mainly intended for Jilian Xiong. The goal is to allow another user besides PM to maintain and quickly troublshoot the daily forecast.

## hyak things

**Routine hyak mainenance** happens every second Tuesday of the month. They shut down hyak from about 9 AM to 5 PM. Any job with for which the SBATCH "time" property would put it past 9 AM will be halted. We use a default of #SBATCH --time=02:00:00 in our sbatch scripts so a job that is started after 7 AM will likely be halted, even if the job only runs an hour.

**You need to check on the klone crontab after each maintenance.** Lately the process has been deleting my crontab and it needs to be reinstalled by hand. The crontab is located in LO/driver/crontabs. This will be a problem for any user other than PM, because I assume they will not be able to set the crontab in my klone account. If this becomes an issue the Operator (JX) may have to set up the forecast to run from her own account on klone. This is something we should practice in advance.

David Darr has given JX **write permission** in /dat1/parker on apogee, where all the pre-and post-processing steps happen.

## post-processing things

**The three post-processing jobs** will run for 24 hours each, waiting patiently until they find all the history files they expect for a given forecast. Thus even if the forecasts don't happen until late in the day, you should not need to do anything about the post-processing.

PM's account on apogee has built in **permission (ssh keys)** allowing it to write files to the APL server, UW kopah storage, and the UW homer computer where the LiveOcean website is. These permissions may require work to replicate for another user. As long as the pre- and post-processing jobs run by cron on apogee in PM's account are able to find the files they need they should work fine.

## HYCOM issues

With **ocn03** if it does not work then it should go to planB, but there can always be new errors that are not yet accounted for, so we have to keep an eye on it. If it persists in being slow I would (i) look for related messages on the HYCOM Forum, and (ii) possibly put a message on the Forum describing the problem.

Also **ocn03** might run so slowly that it goes past 3 AM when the cas7 forecast is supposed to start. If ocn03 eventually works then the 8 AM backup forecast will take care of things automatically. If ocn03 does go to planB it is good practice to try to re-run it by hand later in the day, but not between 9 and 10 AM Pacific time because then the hycom server is being serviced.

It would be good practice for you (Jilian) to **try running ocn03** for the forecast yourself to see if there are any problems. If you did need to do that while I am away then I assume you would do it in your apogee account and then copy the results by hand to my apogee account so that klone will find them.

Also, it would be good in the near future to develop "**ocn04**" which would use GLORYS to get the backfill or forecast forcing and put it in LO format. This would give us a robust backup in the possibility that hycom completely failed.

Similarly we should develop **new atm forcing code** that uses one of the national or international weather models as a backup to our UW WRF output, which may go away someday.

## post-processing problems

Sometimes it can happen that when you get the email about post1.log from apogee it only shows the first job (currently the wgh nest forcing) as having completed, and then nothing else, even though it had plenty of time.

Usually the problem is in line 90 of LO/post/surface1/post_main.py (the next post-processing job), where it is using scp to copy the surface.nc file to the APL server:
```
post_argfun.copy_to_server(Ldir, out_fn)
```
What I have seen in the past is that this will just keep trying, even though the APL server is down.

Things to do:

1. see if you can look at files on the APL server: https://liveocean.apl.uw.edu/output/

2. if the server is indeed having a problem, alert Alex Dioso and Troy Tanner who are in charge of it (adioso@uw.edu, troyt@apl.washington.edu). They will reboot the server and this typically fixes the problem.

3. Then what usually happens is that the post processing will proceed on its own, and no other intervention will be required. The oly1 nest forcing will be created too late for the first forecast, but is usually there for the backup forecast.