// Contre-épreuve : compilé avec les mêmes options, ce code doit contenir une fma.
// Sinon, le test ne prouverait rien (options ou expression régulière inadaptées).

extern "C" double control(double a, double b, double c) { return a * b + c; }
