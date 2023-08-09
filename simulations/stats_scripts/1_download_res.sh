#! /usr/bin/env bash
#
#
# ssh lambda  '
# 	cd /local/chib/toconnor_grp/bing/posseleff_simulations/run4/;
# 	fd -I "\.ne|.member.pq" resdir >tmp.txt
# 	tar -hczf tmp.tgz -T tmp.txt'
# rsync -a lambda:/local/chib/toconnor_grp/bing/posseleff_simulations/run4/tmp.tgz ./
# tar -xf tmp.tgz
# rm -rf tmp.tgz

ssh lambda  '
	cd /local/chib/toconnor_grp/bing/posseleff_simulations/run4/;
	fd -I -p ".*10003.*ibdobj.*|.*10000.*ibdobj.*|.*20003.*ibdobj.*|.*20000.*ibdobj.*" resdir/ > tmp.txt
	tar -hczf tmp.tgz -T tmp.txt'
rsync -a lambda:/local/chib/toconnor_grp/bing/posseleff_simulations/run4/tmp.tgz ./
tar -xf tmp.tgz
rm -rf tmp.tgz
