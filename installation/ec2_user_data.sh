#!/bin/bash
# EC2 user data (paste into the launch template, Advanced details -> User data) that builds
# and tests a cartloader image unattended as soon as the instance boots. Requires the
# template's "Shutdown behavior" to be "Terminate".
#
# Progress and the final report appear in the instance's system log (console: Actions ->
# Monitor and troubleshoot -> Get system log) and in /var/log/cloud-init-output.log.
# Nothing is pushed; to publish, SSH in, `docker login`, and run the push commands at the
# end of the report. Then `sudo shutdown -h now` terminates the instance.

# Safety net: terminate after 24 h even if nobody comes back (cancel: sudo shutdown -c)
shutdown -h +1440

dnf install -y docker git wget unzip tmux
usermod -a -G docker ec2-user
systemctl enable --now docker

# To build another branch or tag, add `git checkout <ref> &&` after `cd cartloader`.
sudo -u ec2-user -i bash -c '
	git clone https://github.com/seqscope/cartloader.git &&
	cd cartloader &&
	bash installation/docker_release.sh --version "$(date +%Y%m%d)a"
'
