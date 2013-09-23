raytra:
	g++ -O2 -I. -I/usr/include/OpenEXR -L/usr/local/lib -lIlmImf -lImath -lHalf -Wall *.cpp -o raytra
