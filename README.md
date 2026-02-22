# Wild Goppa Code Decorder by Peterson Algorithm over Prime Fields
gcc -O3 fujiyama.c  
ulimit -s unlimited  
./a.out  

# 素体上のGoppa符号の利点

・符号長の自由度（GF(2)だと拡大体の大きさに差がありすぎる、素体だと素数のサイズが細かく選べるので符号長は柔軟）

・実装の容易さ（素数で割るだけなので、下手に拡大体で計算するより楽で早い）
