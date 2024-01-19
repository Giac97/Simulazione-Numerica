cp input.circle input.in
./Exercise9
cp best_paths.txt paths_circle.txt
cp distance.out distance_ci.txt

cp input.square input.in
./Exercise9
cp best_paths.txt paths_square.txt
cp distance.out distance_sq.txt

cp input.usa input.in
./Exercise9
cp best_paths.txt paths_usa.txt
cp distance.out distance_us.txt

git add paths_*
git add distance_*
git commit -m "performed circ, sq and usa"
git push
