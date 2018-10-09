# Instruction for setting up regular forecast (launchd in Mac OS)

* Modify the `.plist` file (describe the job) as necessary
  * label: name of the job within launchd
  * Program: the full path to the executable
  * RunAtLoad: run as the job definition is loaded or not
  * StartCalendarInterval: specific time and date

* Store the `.plist` file to relevant directory
  * `~/Library/LaunchAgents` runs on behalf of the logged-in user
  * `/Library/LaunchDaemons` runs on behalf of the root users

* `launchctl list` to see current jobs

* `launchctl load /Library/LaunchDaemons/<*.plist>` to load the job

* `launchctl list | grep local.forecast_OCSPapa` to check the status of the job,

   which returns `<pid> <status> local.forecast_OCSPapa`

   An exit code (<status>) of 0 means either that the program was finished
   successfully or that is has not been started yet. Positive numbers are
   returned for program errors, while negative ones mean that the job has
   terminated due to a signal.

* `launchctl unload ~/Library/LaunchAgents/<*.plist>` to unload the job
