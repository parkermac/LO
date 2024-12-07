# README for Running the LiveOcean Forecast Model

## These notes were started 2024.12.07 by Parker MacCready, and are mainly intended for Jilian Xiong. The goal is to allow another user besides PM to maintain and quickly troublshoot the daily forecast.

Routine hyak mainenance happens every second Tuesday of the month. They shut down hyak from about 9 AM to 5 PM. Any job with for which the SBATCH "time" property would put it past 9 AM will be halted. We use a default of #SBATCH --time=02:00:00 in our sbatch scripts so a job that is started after 7 AM will likely be halted, even if the job only runs an hour.

You need to check on the klone crontab after each maintenance. Lately the process has been deleting my crontab and it needs to be reinstalled by hand. The crontab is located in LO/driver/crontabs. This will be a problem for any user other than PM, because I assume they will not be able to set the crontab in my klone account. If this becomes an issue the Operator (JX) may have to set up the forecast to run from her own account on klone. This is something we should practice in advance.

David Darr has given JX write permission in /dat1/parker on apogee, where all the pre-and post-processing steps happen.

PM's account on apogee has built in permission (ssh keys) allowing it to write files to the APL server, UW kopah storage, and the UW homer computer where the LiveOcean website is. These permissions may require work to replicate for another user. As long as the pre- and post-processing jobs run by cron on apogee in PM's account are able to find the files they need they should work fine.