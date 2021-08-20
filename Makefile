
install: rmath.m
	mmc --make  --c-include-directory `pkg-config --cflags libRmath` `pkg-config --libs libRmath` --no-libgrade --libgrade asm_fast.gc --install-prefix . librmath.install

test: install test_rmath.m
	mmc --make --mld ./lib/mercury --ml rmath `pkg-config --libs libRmath` -E test_rmath && ./test_rmath

PROXY += test install
