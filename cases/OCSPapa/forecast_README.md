# Instruction for setting up regular forecast (via launchd in Mac OS)

* Modify the '.plist' file (describe the job) as necessary
  * label: name of the job within launchd
  * Program: the full path to the executable
  * RunAtLoad: run as the job definition is loaded or not
  * StartCalendarInterval: specific time and date

* Store the '.plist' file to relevant directory
  * “~/Library/LaunchAgents” runs on behalf of the logged-in user
  * “/Library/LaunchDaemons” runs on behalf of the root users

* 'launchctl list' to see current jobs
* 'launchctl load ~/Library/LaunchAgents/<*.plist>' to load the job
* 'launchctl unload ~/Library/LaunchAgents/<*.plist>' to unload the job
