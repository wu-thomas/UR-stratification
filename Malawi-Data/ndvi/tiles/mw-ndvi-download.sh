#!/bin/bash

GREP_OPTIONS=''

cookiejar=$(mktemp cookies.XXXXXXXXXX)
netrc=$(mktemp netrc.XXXXXXXXXX)
chmod 0600 "$cookiejar" "$netrc"
function finish {
  rm -rf "$cookiejar" "$netrc"
}

trap finish EXIT
WGETRC="$wgetrc"

prompt_credentials() {
    echo "Enter your Earthdata Login or other provider supplied credentials"
    read -p "Username (0425tao): " username
    username=${username:-0425tao}
    read -s -p "Password: " password
    echo "machine urs.earthdata.nasa.gov login $username password $password" >> $netrc
    echo
}

exit_with_error() {
    echo
    echo "Unable to Retrieve Data"
    echo
    echo $1
    echo
    echo "https://e4ftl01.cr.usgs.gov//MODV6_Cmp_B/MOLT/MOD13A3.006/2008.12.01/MOD13A3.A2008336.h21v09.006.2015186060356.hdf"
    echo
    exit 1
}

prompt_credentials
  detect_app_approval() {
    approved=`curl -s -b "$cookiejar" -c "$cookiejar" -L --max-redirs 5 --netrc-file "$netrc" https://e4ftl01.cr.usgs.gov//MODV6_Cmp_B/MOLT/MOD13A3.006/2008.12.01/MOD13A3.A2008336.h21v09.006.2015186060356.hdf -w %{http_code} | tail  -1`
    if [ "$approved" -ne "302" ]; then
        # User didn't approve the app. Direct users to approve the app in URS
        exit_with_error "Please ensure that you have authorized the remote application by visiting the link below "
    fi
}

setup_auth_curl() {
    # Firstly, check if it require URS authentication
    status=$(curl -s -z "$(date)" -w %{http_code} https://e4ftl01.cr.usgs.gov//MODV6_Cmp_B/MOLT/MOD13A3.006/2008.12.01/MOD13A3.A2008336.h21v09.006.2015186060356.hdf | tail -1)
    if [[ "$status" -ne "200" && "$status" -ne "304" ]]; then
        # URS authentication is required. Now further check if the application/remote service is approved.
        detect_app_approval
    fi
}

setup_auth_wget() {
    # The safest way to auth via curl is netrc. Note: there's no checking or feedback
    # if login is unsuccessful
    touch ~/.netrc
    chmod 0600 ~/.netrc
    credentials=$(grep 'machine urs.earthdata.nasa.gov' ~/.netrc)
    if [ -z "$credentials" ]; then
        cat "$netrc" >> ~/.netrc
    fi
}

fetch_urls() {
  if command -v curl >/dev/null 2>&1; then
      setup_auth_curl
      while read -r line; do
        # Get everything after the last '/'
        filename="${line##*/}"

        # Strip everything after '?'
        stripped_query_params="${filename%%\?*}"

        curl -f -b "$cookiejar" -c "$cookiejar" -L --netrc-file "$netrc" -g -o $stripped_query_params -- $line && echo || exit_with_error "Command failed with error. Please retrieve the data manually."
      done;
  elif command -v wget >/dev/null 2>&1; then
      # We can't use wget to poke provider server to get info whether or not URS was integrated without download at least one of the files.
      echo
      echo "WARNING: Can't find curl, use wget instead."
      echo "WARNING: Script may not correctly identify Earthdata Login integrations."
      echo
      setup_auth_wget
      while read -r line; do
        # Get everything after the last '/'
        filename="${line##*/}"

        # Strip everything after '?'
        stripped_query_params="${filename%%\?*}"

        wget --load-cookies "$cookiejar" --save-cookies "$cookiejar" --output-document $stripped_query_params --keep-session-cookies -- $line && echo || exit_with_error "Command failed with error. Please retrieve the data manually."
      done;
  else
      exit_with_error "Error: Could not find a command-line downloader.  Please install curl or wget"
  fi
}

fetch_urls <<'EDSCEOF'
https://e4ftl01.cr.usgs.gov//MODV6_Cmp_B/MOLT/MOD13A3.006/2008.12.01/MOD13A3.A2008336.h21v09.006.2015186060356.hdf
https://e4ftl01.cr.usgs.gov//MODV6_Cmp_B/MOLT/MOD13A3.006/2008.12.01/MOD13A3.A2008336.h21v10.006.2015186060413.hdf
https://e4ftl01.cr.usgs.gov//MODV6_Cmp_B/MOLT/MOD13A3.006/2008.11.01/MOD13A3.A2008306.h21v10.006.2015181145917.hdf
https://e4ftl01.cr.usgs.gov//MODV6_Cmp_B/MOLT/MOD13A3.006/2008.11.01/MOD13A3.A2008306.h21v09.006.2015181145923.hdf
https://e4ftl01.cr.usgs.gov//MODV6_Cmp_B/MOLT/MOD13A3.006/2008.10.01/MOD13A3.A2008275.h21v10.006.2015181145106.hdf
https://e4ftl01.cr.usgs.gov//MODV6_Cmp_B/MOLT/MOD13A3.006/2008.10.01/MOD13A3.A2008275.h21v09.006.2015181144221.hdf
https://e4ftl01.cr.usgs.gov//MODV6_Cmp_B/MOLT/MOD13A3.006/2008.09.01/MOD13A3.A2008245.h21v09.006.2015180222406.hdf
https://e4ftl01.cr.usgs.gov//MODV6_Cmp_B/MOLT/MOD13A3.006/2008.09.01/MOD13A3.A2008245.h21v10.006.2015180222356.hdf
https://e4ftl01.cr.usgs.gov//MODV6_Cmp_B/MOLT/MOD13A3.006/2008.08.01/MOD13A3.A2008214.h21v09.006.2015179232108.hdf
https://e4ftl01.cr.usgs.gov//MODV6_Cmp_B/MOLT/MOD13A3.006/2008.08.01/MOD13A3.A2008214.h21v10.006.2015179232104.hdf
https://e4ftl01.cr.usgs.gov//MODV6_Cmp_B/MOLT/MOD13A3.006/2008.07.01/MOD13A3.A2008183.h21v10.006.2015177045508.hdf
https://e4ftl01.cr.usgs.gov//MODV6_Cmp_B/MOLT/MOD13A3.006/2008.07.01/MOD13A3.A2008183.h21v09.006.2015177045511.hdf
https://e4ftl01.cr.usgs.gov//MODV6_Cmp_B/MOLT/MOD13A3.006/2008.06.01/MOD13A3.A2008153.h21v10.006.2015176071640.hdf
https://e4ftl01.cr.usgs.gov//MODV6_Cmp_B/MOLT/MOD13A3.006/2008.06.01/MOD13A3.A2008153.h21v09.006.2015176072448.hdf
https://e4ftl01.cr.usgs.gov//MODV6_Cmp_B/MOLT/MOD13A3.006/2008.05.01/MOD13A3.A2008122.h21v10.006.2015175110740.hdf
https://e4ftl01.cr.usgs.gov//MODV6_Cmp_B/MOLT/MOD13A3.006/2008.05.01/MOD13A3.A2008122.h21v09.006.2015175114120.hdf
https://e4ftl01.cr.usgs.gov//MODV6_Cmp_B/MOLT/MOD13A3.006/2008.04.01/MOD13A3.A2008092.h21v10.006.2015173161052.hdf
https://e4ftl01.cr.usgs.gov//MODV6_Cmp_B/MOLT/MOD13A3.006/2008.04.01/MOD13A3.A2008092.h21v09.006.2015173161146.hdf
https://e4ftl01.cr.usgs.gov//MODV6_Cmp_B/MOLT/MOD13A3.006/2008.03.01/MOD13A3.A2008061.h21v10.006.2015173023704.hdf
https://e4ftl01.cr.usgs.gov//MODV6_Cmp_B/MOLT/MOD13A3.006/2008.03.01/MOD13A3.A2008061.h21v09.006.2015173025325.hdf
https://e4ftl01.cr.usgs.gov//MODV6_Cmp_B/MOLT/MOD13A3.006/2008.02.01/MOD13A3.A2008032.h21v10.006.2015172131755.hdf
https://e4ftl01.cr.usgs.gov//MODV6_Cmp_B/MOLT/MOD13A3.006/2008.02.01/MOD13A3.A2008032.h21v09.006.2015172131545.hdf
https://e4ftl01.cr.usgs.gov//MODV6_Cmp_B/MOLT/MOD13A3.006/2008.01.01/MOD13A3.A2008001.h21v09.006.2015170035845.hdf
https://e4ftl01.cr.usgs.gov//MODV6_Cmp_B/MOLT/MOD13A3.006/2008.01.01/MOD13A3.A2008001.h21v10.006.2015170042048.hdf
EDSCEOF