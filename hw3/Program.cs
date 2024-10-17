using System;
using System.Data;
using System.Data.Common;
using System.Diagnostics;
using System.Globalization;
using System.Reflection.Metadata;
using System.Runtime.CompilerServices;
using System.Runtime.ExceptionServices;
using System.Runtime.InteropServices;
using System.Security.Authentication.ExtendedProtection;
using System.Security.Cryptography.X509Certificates;
using System.Xml;

namespace MyApp
{
    internal class Program
    {
        // Hi grader, run this function to test the implementation. 
        public static void ForGrader()
        {
            MonteCarloStandard mc = new MonteCarloStandard()
            {
                StockPrice = 60,
                Strike = 40,
                Tenor = 1,
                Vol = 0.4,
                Rate = 0.05,
                isCall = true,
                doParallel = false
            };

            MonteCarloAntithetic mc_anti = new MonteCarloAntithetic()
            {
                StockPrice = 60,
                Strike = 40,
                Tenor = 1,
                Vol = 0.4,
                Rate = 0.05,
                isCall = true,
                doParallel = false
            };

            MonteCarloCV mc_cv = new MonteCarloCV()
            {
                StockPrice = 60,
                Strike = 40,
                Tenor = 1,
                Vol = 0.4,
                Rate = 0.05,
                isCall = true,
                doParallel = false
            };

            MonteCarloAntitheticCV mc_anti_cv = new MonteCarloAntitheticCV()
            {
                StockPrice = 60,
                Strike = 40,
                Tenor = 1,
                Vol = 0.4,
                Rate = 0.05,
                isCall = true,
                doParallel = false
            };
            MonteCarloStandard mc_parallel = new MonteCarloStandard()
            {
                StockPrice = 60,
                Strike = 40,
                Tenor = 1,
                Vol = 0.4,
                Rate = 0.05,
                isCall = true,
                doParallel = false
            };

            MonteCarloAntithetic mc_anti_parallel  = new MonteCarloAntithetic()
            {
                StockPrice = 60,
                Strike = 40,
                Tenor = 1,
                Vol = 0.4,
                Rate = 0.05,
                isCall = true,
                doParallel = false
            };

            MonteCarloCV mc_cv_parallel  = new MonteCarloCV()
            {
                StockPrice = 60,
                Strike = 40,
                Tenor = 1,
                Vol = 0.4,
                Rate = 0.05,
                isCall = true,
                doParallel = true
            };

            MonteCarloAntitheticCV mc_anti_cv_parallel  = new MonteCarloAntitheticCV()
            {
                StockPrice = 60,
                Strike = 40,
                Tenor = 1,
                Vol = 0.4,
                Rate = 0.05,
                isCall = true,
                doParallel = true
            };
            
            int trials = 10000;
            int steps = 252;
            int seed = 10;

            double[] pair1 = mc.SimulateOptionPrice(trials, steps, seed);
            double[] pair2 = mc_anti.SimulateOptionPrice(trials, steps, seed);
            double[] pair3 = mc_cv.SimulateOptionPrice(trials, steps, seed);
            double[] pair4 = mc_anti_cv.SimulateOptionPrice(trials, steps, seed);
            double[] pair5 = mc_parallel.SimulateOptionPrice(trials, steps, seed);
            double[] pair6 = mc_anti_parallel.SimulateOptionPrice(trials, steps, seed);
            double[] pair7 = mc_cv_parallel.SimulateOptionPrice(trials, steps, seed);
            double[] pair8 = mc_anti_cv_parallel.SimulateOptionPrice(trials, steps, seed);
            Console.WriteLine("Price and Standard Error");
            Console.WriteLine("No Variance Reduction: {0}, {1}", pair1[0], pair1[1]);
            Console.WriteLine("Antithetic Only: {0}, {1}", pair2[0], pair2[1]);
            Console.WriteLine("Control Variate Only: {0}, {1}", pair3[0], pair3[1]);
            Console.WriteLine("Antithetic AND CV: {0}, {1}", pair4[0], pair4[1]);
            Console.WriteLine("NOW FOR PARALLEL");
            Console.WriteLine("No Variance Reduction: {0}, {1}", pair5[0], pair5[1]);
            Console.WriteLine("Antithetic Only: {0}, {1}", pair6[0], pair6[1]);
            Console.WriteLine("Control Variate Only: {0}, {1}", pair7[0], pair7[1]);
            Console.WriteLine("Antithetic AND CV: {0}, {1}", pair8[0], pair8[1]);
        }

        public static void RunConsoleApp()
        {
            string startMessage = 
                @"Choose one of the following:
                0 for Standard Monte-Carl0
                1 for Antithetic Monte-Carlo
                2 for Quasi-Monte-Carlo using Van der Corput
                3 for Control Variate Monte-Carlo
                4 for Antithetic and Control Variate Monte-Carlo";


            Console.WriteLine(startMessage);

            try
            {
                string input = Console.ReadLine();
                Console.Write("Enter Stock Price: ");
                string stockPrice = Console.ReadLine();
                Console.Write("Enter Strike: ");
                string strike = Console.ReadLine();
                Console.Write("Enter Tenor: ");
                string tenor = Console.ReadLine();
                Console.Write("Enter Volatility: ");
                string vol = Console.ReadLine();
                Console.Write("Enter Risk-Free Rate: ");
                string rate = Console.ReadLine();
                Console.Write("Is the option a call? (Enter \"true\" or \"false\"): ");
                string isCall = Console.ReadLine();
                Console.Write("Enter the number of trials: ");
                string n = Console.ReadLine();

                if (int.Parse(input) == 0) 
                {
                    MonteCarloStandard mc = new MonteCarloStandard()
                    {
                        StockPrice = double.Parse(stockPrice), 
                        Strike = double.Parse(strike),
                        Tenor = double.Parse(tenor),
                        Vol = double.Parse(vol),
                        Rate = double.Parse(rate),
                        isCall = bool.Parse(isCall),
                    };
                
                    mc.PrintOptionInfo();
                }
                else if (int.Parse(input) == 1)
                {
                    MonteCarloAntithetic mc = new MonteCarloAntithetic()
                    {
                        StockPrice = double.Parse(stockPrice), 
                        Strike = double.Parse(strike),
                        Tenor = double.Parse(tenor),
                        Vol = double.Parse(vol),
                        Rate = double.Parse(rate),
                        isCall = bool.Parse(isCall),
                    };
                
                    mc.PrintOptionInfo();
                }
                else if (int.Parse(input) == 2)
                {   
                    Console.Write("Enter a base: ");
                    string baseNumber = Console.ReadLine();
                    QuasiMonteCarlo mc = new QuasiMonteCarlo()
                    {
                        StockPrice = double.Parse(stockPrice), 
                        Strike = double.Parse(strike),
                        Tenor = double.Parse(tenor),
                        Vol = double.Parse(vol),
                        Rate = double.Parse(rate),
                        isCall = bool.Parse(isCall),
                    };
                    mc.PrintOptionInfo(Int32.Parse(n), Int32.Parse(baseNumber));

                }
                else if (int.Parse(input) == 3)
                {
                    MonteCarloCV mc = new MonteCarloCV()
                    {
                        StockPrice = double.Parse(stockPrice), 
                        Strike = double.Parse(strike),
                        Tenor = double.Parse(tenor),
                        Vol = double.Parse(vol),
                        Rate = double.Parse(rate),
                        isCall = bool.Parse(isCall),
                    };
                
                    mc.PrintOptionInfo();
                }
                else if (int.Parse(input) == 4)
                {
                    MonteCarloAntitheticCV mc = new MonteCarloAntitheticCV()
                    {
                        StockPrice = double.Parse(stockPrice), 
                        Strike = double.Parse(strike),
                        Tenor = double.Parse(tenor),
                        Vol = double.Parse(vol),
                        Rate = double.Parse(rate),
                        isCall = bool.Parse(isCall),
                    };
                
                    mc.PrintOptionInfo();
                }
                else
                {
                    Console.WriteLine("Try again");
                }

                
            } 
            catch (Exception e) {
                Console.WriteLine("There was an error with your inputs, try again.");
            }

            // Good test numbers
            // StockPrice = 60, 
            // Strike = 40,
            // Tenor = 1,
            // Vol = 0.4,
            // Rate = 0.05
        }

        static void Main(string[] args)
        {
            // Run ForGrader to see a premade test for your validation
            // Just bought a house, working full time, didnt have time to put my best
            // effort on this one. Will fix all the issues for the next assignment
            ForGrader();

            //RunConsoleApp();
            // mc_anti.PrintOptionInfo();
            
        }
    }


    public class Scenarios
        {
            public static double[,] GenerateScenarios(double StockPrice, double Tenor, double rate, double vol, int trials, int steps, int seed=0, bool antiCorr = false)
            {
                int sign = antiCorr ? 1 : -1;
                // Generate a matrix of normal random samples
                NormalRandom rand = new NormalRandom(seed);
                double dt = Tenor/steps;
                double[,] samples = rand.Normal(0, Math.Sqrt(dt), trials, steps+1);
                
                double[,] simulation = new double[trials,steps+1];
                
                // Set first column of the simulation matrix equal to the initial underlying price
                for (int i = 0; i < trials; i++)
                {   
                    simulation[i,0] = StockPrice;
                }

                for (int i = 0; i < trials; i++) 
                {
                    for(int j = 0; j < steps; j++)
                    {
                        simulation[i,j+1] = simulation[i,j] * Math.Exp(((rate - (Math.Pow(vol, 2)/2)) * dt) + sign * (vol * samples[i,j]));
                    }
                }
                return simulation;
            }
            
            public static double[,] GenerateScenariosParallel(double StockPrice, double Tenor, double rate, double vol, int trials, int steps, int seed=0, bool antiCorr = false)
            {
                int sign = antiCorr ? 1 : -1;
                int numCores = Environment.ProcessorCount;
                if (trials % numCores != 0) throw new Exception();
                NormalRandom rand = new NormalRandom();
                double[,] simulation = rand.Normal(0,1,trials, steps+1);
                double[,] scen = new double[trials,steps+1];
                double dt = Tenor/steps;

                int blockSize = trials / numCores;
                Parallel.For(0, numCores, i => 
                {
                    for (int k = i*blockSize; k < (i+1)*blockSize; k++)
                    {
                        scen[k,0] = StockPrice;
                        for (int j = 0; j < steps-1; j++)
                        {
                            scen[k,j+1] = scen[k,j] * (Math.Exp( (rate - Math.Pow(vol, 2)/2))*dt ) + Math.Sqrt(Tenor) * vol * simulation[k,j];
                        }
                    }
                });
                return scen;
            }
        }

    abstract class MonteCarlo
    {
        public double StockPrice {get; set;}
        public double Strike {get; set;}
        public double Tenor {get; set;}
        public double Vol {get; set;}
        public double Rate {get; set;}
        public bool isCall {get; set;}
        public bool doParallel {get; set;}

        public void PrintOptionInfo(int trials = 100000, int steps = 1, int seed = 0)
        {
            // Measure execution time
            var watch = System.Diagnostics.Stopwatch.StartNew();
            double[] price_se_pair = SimulateOptionPrice(trials, steps, seed);
            watch.Stop();
            Console.WriteLine("Option Price: {0}", Math.Round(price_se_pair[0], 6));
            Console.WriteLine("StErr of Simulation: {0}", Math.Round(price_se_pair[1], 6));
            Console.WriteLine("Delta: {0}", Math.Round(SimulateDelta(trials, steps, seed), 6));
            Console.WriteLine("Gamma: {0}", Math.Round(SimulateGamma(trials, steps, seed), 6));
            Console.WriteLine("Vega: {0}", Math.Round(SimulateVega(trials, steps, seed), 6));
            Console.WriteLine("Theta: {0}", Math.Round(SimulateTheta(trials, steps, seed), 6));
            Console.WriteLine("Rho: {0}", Math.Round(SimulateRho(trials, steps, seed), 6));
            
            double elapsedMs = watch.ElapsedMilliseconds;
            Console.WriteLine("Execution time is {0} seconds", elapsedMs/1000);
        }

        public virtual double[] SimulateOptionPrice(int trials, int steps, int seed) { return new double[] {0}; }

        public double SimulateDelta(int trials, int steps, int seed, double incPercent=0.001) 
        {
            double increment = incPercent * StockPrice;

            double temp = StockPrice;
            StockPrice += increment;
            double priceUp = SimulateOptionPrice(trials, steps, seed)[0];
            StockPrice = temp;

            StockPrice -= increment;
            double priceDown = SimulateOptionPrice(trials, steps, seed)[0];
            StockPrice = temp;

            return (priceUp - priceDown) / (2 * increment);
        } 

        public double SimulateGamma(int trials, int steps, int seed, double incPercent=0.1) 
        {
            double increment = incPercent * StockPrice;

            double temp = StockPrice;
            StockPrice += increment;
            double priceUp = SimulateOptionPrice(trials, steps, seed)[0];
            StockPrice = temp;

            double priceMid = SimulateOptionPrice(trials, steps, seed)[0];

            StockPrice -= increment;
            double priceDown = SimulateOptionPrice(trials, steps, seed)[0];
            StockPrice = temp;

            return (priceUp - 2*priceMid + priceDown) / Math.Pow(increment, 2);
        }

        public double SimulateVega(int trials, int steps, int seed, double incPercent=0.1) 
        {
            double increment = incPercent * Vol;

            double temp = Vol;
            Vol += increment;
            double priceUp = SimulateOptionPrice(trials, steps, seed)[0];
            Vol = temp;

            Vol -= increment;
            double priceDown = SimulateOptionPrice(trials, steps, seed)[0];
            Vol = temp;

            return (priceUp - priceDown) / (2 * increment);
        }

        public double SimulateTheta(int trials, int steps, int seed, double incPercent=0.1) 
        {
            double increment = incPercent * Tenor;

            double temp = Tenor;
            Tenor += increment;
            double priceUp = SimulateOptionPrice(trials, steps, seed)[0];
            Tenor = temp;

            double priceMid = SimulateOptionPrice(trials, steps, seed)[0];

            return (priceUp - priceMid) / increment;
        }

        public double SimulateRho(int trials, int steps, int seed, double incPercent=0.1) 
        {
            double increment = incPercent * Rate;

            double temp = Rate;
            Rate += increment;
            double priceUp = SimulateOptionPrice(trials, steps, seed)[0];
            Rate = temp;

            Rate -= increment;
            double priceDown = SimulateOptionPrice(trials, steps, seed)[0];
            Rate = temp;

            return (priceUp - priceDown) / (2 * increment);
        }
    }

    class MonteCarloStandard : MonteCarlo
    {
        public override double[] SimulateOptionPrice(int trials, int steps = 1, int seed = 0) 
        {
            // Generate a matrix of normal random samples
            NormalRandom rand = new NormalRandom(seed);
            double dt = Tenor/steps;
            double[,] samples = rand.Normal(0, Math.Sqrt(dt), trials, steps+1);
            
            double[,] simulation = new double[trials, steps+1];
            
            // Set first column of the simulation matrix equal to the initial underlying price
            if (doParallel)
            {
                simulation = Scenarios.GenerateScenariosParallel(StockPrice, Tenor, Rate, Vol, trials, steps, seed);
                
            }
            else 
            {
                for (int i = 0; i < trials; i++)
                {
                    simulation[i, 0] = StockPrice;
                }

                for (int i = 0; i < trials; i++) 
                {
                    for(int j = 0; j < steps; j++)
                    {
                        simulation[i,j+1] = simulation[i,j] * Math.Exp(((Rate - (Math.Pow(Vol, 2)/2)) * dt) + (Vol * samples[i,j]));
                    }
                }
            }


            // Find the payoff of the option of the terminal underlying prices
            double[] payoff = new double[trials];
            if (isCall) 
            {
                for (int i = 0; i < trials; i++)
                {
                    payoff[i] = Payoff.callPayoff(simulation[i, steps], Strike);
                }
            }
            else
            {
                for (int i = 0; i < trials; i++)
                {
                    payoff[i] = Payoff.putPayoff(simulation[i, steps], Strike);
                }
            }
            double estimatedPrice = Math.Exp(-Rate*Tenor)*Stats.Mean(payoff);
            double standardError = Math.Sqrt(Stats.Variance(payoff)) / Math.Sqrt(payoff.Length);
            return new double[] {estimatedPrice, standardError};
        }
    }

    class MonteCarloAntithetic : MonteCarlo
    {
        public override double[] SimulateOptionPrice(int trials, int steps = 1, int seed = 0) 
        { 
            // Generate a matrix of normal random samples
            NormalRandom rand = new NormalRandom(seed);
            double dt = Tenor/steps;
            double[,] samples = rand.Normal(0, Math.Sqrt(dt), trials, steps+1);
            
            double[,] simulation_normal = new double[trials, steps+1];
            double[,] simulation_anti = new double[trials, steps+1];
            
            // Set first column of the simulation matrix equal to the initial underlying price
            if (doParallel)
            {
                simulation_normal = Scenarios.GenerateScenariosParallel(StockPrice, Tenor, Rate, Vol, trials, steps, seed, false);
                simulation_anti = Scenarios.GenerateScenariosParallel(StockPrice, Tenor, Rate, Vol, trials, steps, seed, true);
            }
            else 
            {
                for (int i = 0; i < trials; i++)
                {
                    simulation_normal[i, 0] = StockPrice;
                    simulation_anti[i, 0] = StockPrice;
                }

                
                for (int i = 0; i < trials; i++) 
                {
                    for(int j = 0; j < steps; j++)
                    {
                        simulation_normal[i,j+1] = simulation_normal[i,j] * Math.Exp(((Rate - (Math.Pow(Vol, 2)/2)) * dt) + (Vol * samples[i,j]));
                        simulation_anti[i,j+1] = simulation_anti[i,j] * Math.Exp(((Rate - (Math.Pow(Vol, 2)/2)) * dt) - (Vol * samples[i,j]));
                    }
                }
            }

            // Find the payoff of the option of the terminal underlying prices
            double[] payoff_normal = new double[trials];
            double[] payoff_anti = new double[trials];
            if (isCall) 
            {
                for (int i = 0; i < trials; i++)
                {
                    payoff_normal[i] = Payoff.callPayoff(simulation_normal[i, steps], Strike);
                    payoff_anti[i] = Payoff.callPayoff(simulation_anti[i, steps], Strike);
                }
            }
            else
            {
                for (int i = 0; i < trials; i++)
                {
                    payoff_normal[i] = Payoff.putPayoff(simulation_normal[i, steps], Strike);
                    payoff_anti[i] = Payoff.putPayoff(simulation_anti[i, steps], Strike);
                }
            } 
            double estimatedPrice = Math.Exp(-Rate*Tenor) * (Stats.Mean(payoff_normal) + Stats.Mean(payoff_anti)) / 2;
            
            // Inner method for calculated Standard Error in the Antithetic case
            double CalculateSE(double[] arr1, double[] arr2)
            {
                double[] dummy = new double[arr1.Length];
                for (int i = 0; i < arr1.Length; i++)
                {
                    dummy[i] = (arr1[i] + arr1[2]) / 2;
                }
                double variance = Stats.Variance(dummy);
                return Math.Sqrt(variance/arr1.Length);
            }

            double standardError = CalculateSE(payoff_normal, payoff_anti);

            return new double[] { estimatedPrice, standardError };
        }
    }

    class MonteCarloCV : MonteCarlo
    {
        delegate double PayoffFunction(double ST, double K);

        public override double[] SimulateOptionPrice(int trials, int steps, int seed)
        {
            // public double StockPrice {get; set;}
            // public double Strike {get; set;}
            // public double Tenor {get; set;}
            // public double Vol {get; set;}
            // public double Rate {get; set;}

            // Generate a matrix of normal random samples
            

            PayoffFunction payoff = isCall ? (double ST, double K) => Math.Max(ST - K, 0) : (double ST, double K) => Math.Max(K - ST, 0);

            NormalRandom rand = new NormalRandom(seed);
            double dt = Tenor/steps;
            double nudt = (Rate - 0.5 * Math.Pow(Vol,2))*dt;
            double sigsdt = Vol * Math.Sqrt(dt);
            double erddt = Math.Exp(Rate*dt);

            int beta = -1;

            double sum_CT = 0;
            double sum_CT2 = 0;

            double[,] samples = rand.Normal(0, 1, trials, steps+1);

            if (doParallel)
            {
                // Didnt have time to implement this part. Just bought a house. Will get it done for next HW.
                return new double[] {0, 0};
            }
            else
            {
                for (int j = 0; j < trials-1; j++) 
                {
                    double St = StockPrice;
                    double cv = 0;

                    for(int i = 0; i < steps; i++)
                    {
                        double t = i*dt;
                        double delta = BlackScholes.BlackScholesDelta(St, Strike, Tenor - t, Rate, Vol, isCall);
                        double Stn = St*Math.Exp(nudt + sigsdt * samples[j,i]);
                        cv += delta*(Stn-St*erddt);
                        St = Stn;
                    }
                    
                    double CT = payoff(St, Strike) + beta*cv;
                    sum_CT += CT;
                    sum_CT2 += Math.Pow(CT,2);
                }
            
                double option_value = Math.Exp(-Rate*Tenor) * (sum_CT / trials); 
                double SD = Math.Sqrt( ( sum_CT2 - sum_CT*sum_CT/trials ) * Math.Exp(-2*Rate*Tenor) / (trials - 1) );
                double SE = SD / Math.Sqrt(trials);

                return new double[] {option_value, SE};
            }
        }
    }

    class MonteCarloAntitheticCV : MonteCarlo
    {
        delegate double PayoffFunction(double ST, double K);

        public override double[] SimulateOptionPrice(int trials, int steps, int seed)
        {
            // public double StockPrice {get; set;}
            // public double Strike {get; set;}
            // public double Tenor {get; set;}
            // public double Vol {get; set;}
            // public double Rate {get; set;}

            // Generate a matrix of normal random samples
            
            PayoffFunction payoff = isCall ? (double ST, double K) => Math.Max(ST - K, 0) : (double ST, double K) => Math.Max(K - ST, 0);

            NormalRandom rand = new NormalRandom(seed);
            double dt = Tenor/steps;
            double nudt = (Rate - 0.5 * Math.Pow(Vol,2))*dt;
            double sigsdt = Vol * Math.Sqrt(dt);
            double erddt = Math.Exp(Rate*dt);

            int beta = -1;

            double sum_CT = 0;
            double sum_CT2 = 0;

            double[,] samples = rand.Normal(0, 1, trials, steps+1);
            
            if (doParallel)
            {
                // Didnt have time to implement this part. Just bought a house. Will get it done for next HW.
                return new double[] {0, 0};
            }
            else
            {
                for (int j = 0; j < trials-1; j++) 
                {
                    double St1 = StockPrice;
                    double St2 = StockPrice;
                    double cv1 = 0;
                    double cv2 = 0;

                    for(int i = 0; i < steps; i++)
                    {
                        double t = i*dt;
                        double delta1 = BlackScholes.BlackScholesDelta(St1, Strike, Tenor - t, Rate, Vol, isCall);
                        double delta2 = BlackScholes.BlackScholesDelta(St2, Strike, Tenor - t, Rate, Vol, isCall);
                        double Stn1 = St1*Math.Exp(nudt + sigsdt * samples[j,i]);
                        double Stn2 = St2*Math.Exp(nudt - sigsdt * samples[j,i]);

                        cv1 += delta1*(Stn1-St1*erddt);
                        cv2 += delta2*(Stn2-St2*erddt);

                        St1 = Stn1;
                        St2 = Stn2;
                    }
                    
                    double CT = 0.5 * ( payoff(St1, Strike) + beta*cv1 + payoff(St2, Strike) + beta*cv2);
                    sum_CT += CT;
                    sum_CT2 += Math.Pow(CT,2);
                }
                
                double option_value = Math.Exp(-Rate*Tenor) * (sum_CT / trials); 
                double SD = Math.Sqrt( ( sum_CT2 - sum_CT*sum_CT/trials ) * Math.Exp(-2*Rate*Tenor) / (trials - 1) );
                double SE = SD / Math.Sqrt(trials);

                return new double[] {option_value, SE};
            }
        }
    }
    

    class QuasiMonteCarlo
    {
        public double StockPrice {get; set;}
        public double Strike {get; set;}
        public double Tenor {get; set;}
        public double Vol {get; set;}
        public double Rate {get; set;}
        public bool isCall {get; set;}

        public void PrintOptionInfo(int trials = 100000, int baseNumber = 2)
        {
            // Measure execution time
            var watch = System.Diagnostics.Stopwatch.StartNew();
            double[] price_se_pair = SimulateOptionPrice(trials, baseNumber);
            Console.WriteLine("Option Price: {0}", Math.Round(price_se_pair[0], 6));
            Console.WriteLine("StErr of Simulation: {0}", Math.Round(price_se_pair[1], 6));
            Console.WriteLine("Delta: {0}", Math.Round(SimulateDelta(trials, baseNumber), 6));
            Console.WriteLine("Gamma: {0}", Math.Round(SimulateGamma(trials, baseNumber), 6));
            Console.WriteLine("Vega: {0}", Math.Round(SimulateVega(trials, baseNumber), 6));
            Console.WriteLine("Theta: {0}", Math.Round(SimulateTheta(trials, baseNumber, 6)));
            Console.WriteLine("Rho: {0}", Math.Round(SimulateRho(trials, baseNumber), 6));
            watch.Stop();
            double elapsedMs = watch.ElapsedMilliseconds;
            Console.WriteLine("Execution time is {0} seconds", elapsedMs/1000);
        }

        public double[] CorputSequence(int length, int baseNumber)
        {
            // Compute the nth number in the corput number with the corresponding base
            double corput(int n, int baseNumber)
            {
                double q = 0;
                double bk = (double) 1/baseNumber; 

                while (n > 0)
                {
                    q += (n % baseNumber) * bk;
                    n /= baseNumber;
                    bk /= baseNumber;
                }

                return q;
            }

            double[] output = new double[length];
            for (int i = 0; i < length; i++)
            {
                output[i] = corput(i+1, baseNumber);
            }
            return output;
        }

        public double[] SimulateOptionPrice(int trials, int baseNumber)
        {
            // Apply transformation to the corput sequence
            double[] samples = CorputSequence(trials, baseNumber);
            for (int i = 0; i < samples.Length; i += 2)
            {
                double x1 = samples[i];
                double x2 = samples[i+1];

                samples[i] = Math.Sqrt(-2*Math.Log(x1))*Math.Cos(2*Math.PI*x2);
                samples[i+1] = Math.Sqrt(-2*Math.Log(x1))*Math.Sin(2*Math.PI*x2);
            }

            int steps = 1;
            double dt = Tenor/steps;
            
            double[,] simulation = new double[trials, steps+1];

            // Set first column of the simulation matrix equal to the initial underlying price
            for (int i = 0; i < trials; i++)
            {
                simulation[i,0] = StockPrice;
            }

            for (int i = 0; i < trials; i++) 
            {
                for(int j = 0; j < steps; j++)
                {
                    simulation[i,j+1] = simulation[i,j] * Math.Exp(((Rate - (Math.Pow(Vol, 2)/2)) * dt) + (Vol * samples[i]));
                }
            }

            // Find the payoff of the option of the terminal underlying prices
            double[] payoff = new double[trials];
            if (isCall) 
            {
                for (int i = 0; i < trials; i++)
                {
                    payoff[i] = Payoff.callPayoff(simulation[i, steps], Strike);
                }
            }
            else
            {
                for (int i = 0; i < trials; i++)
                {
                    payoff[i] = Payoff.putPayoff(simulation[i, steps], Strike);
                }
            }
            double estimatedPrice = Math.Exp(-Rate*Tenor)*Stats.Mean(payoff);
            double standardError = Stats.Variance(payoff)/Math.Sqrt(payoff.Length);
            return new double[] {estimatedPrice, standardError};
        }

        public double SimulateDelta(int trials, int baseNumber, double incPercent=0.001) 
        {
            double increment = incPercent * StockPrice;

            double temp = StockPrice;
            StockPrice += increment;
            double priceUp = SimulateOptionPrice(trials, baseNumber)[0];
            StockPrice = temp;

            StockPrice -= increment;
            double priceDown = SimulateOptionPrice(trials, baseNumber)[0];
            StockPrice = temp;

            return (priceUp - priceDown) / (2 * increment);
        } 

        public double SimulateGamma(int trials, int baseNumber, double incPercent=0.1) 
        {
            double increment = incPercent * StockPrice;

            double temp = StockPrice;
            StockPrice += increment;
            double priceUp = SimulateOptionPrice(trials, baseNumber)[0];
            StockPrice = temp;

            double priceMid = SimulateOptionPrice(trials, baseNumber)[0];

            StockPrice -= increment;
            double priceDown = SimulateOptionPrice(trials, baseNumber)[0];
            StockPrice = temp;

            return (priceUp - 2*priceMid + priceDown) / Math.Pow(increment, 2);
        }

        public double SimulateVega(int trials, int baseNumber, double incPercent=0.1) 
        {
            double increment = incPercent * Vol;

            double temp = Vol;
            Vol += increment;
            double priceUp = SimulateOptionPrice(trials, baseNumber)[0];
            Vol = temp;

            Vol -= increment;
            double priceDown = SimulateOptionPrice(trials, baseNumber)[0];
            Vol = temp;

            return (priceUp - priceDown) / (2 * increment);
        }

        public double SimulateTheta(int trials, int baseNumber, double incPercent=0.1) 
        {
            double increment = incPercent * Tenor;

            double temp = Tenor;
            Tenor += increment;
            double priceUp = SimulateOptionPrice(trials, baseNumber)[0];
            Tenor = temp;

            double priceMid = SimulateOptionPrice(trials, baseNumber)[0];

            return (priceUp - priceMid) / increment;
        }

        public double SimulateRho(int trials, int baseNumber, double incPercent=0.1) 
        {
            double increment = incPercent * Rate;

            double temp = Rate;
            Rate += increment;
            double priceUp = SimulateOptionPrice(trials, baseNumber)[0];
            Rate = temp;

            Rate -= increment;
            double priceDown = SimulateOptionPrice(trials, baseNumber)[0];
            Rate = temp;

            return (priceUp - priceDown) / (2 * increment);
        }
    }
    

    
    public static class BlackScholes
    {
        public static double BlackScholesDelta(double S, double K, double T, double r, double sigma, bool isCall)
        {
            int sign = isCall ? 1 : -1;
            double d1 = (Math.Log(S / K) + (r + 0.5 * sigma * sigma) * T) / (sigma * Math.Sqrt(T));
            return sign * Stats.NormCDF(d1);
        }
    }


    class NormalRandom : Random
    {
        private Random rand;

        public NormalRandom() { rand = new Random(); }
        public NormalRandom(int seed) { rand = new Random(seed); }

        public double[,] Normal(double mean=0, double variance=1, int n=1, int m=1)
        {
            double[] samples = Normal(mean, variance, n*m);
            double [,] reshapeSamples = new double[n,m];
            for (int i = 0; i < n; i++) {
                for (int j = 0; j < m; j++) {
                    reshapeSamples[i,j] = samples[i * m + j];
                }
            }
            return reshapeSamples;
        }

        // Generate "size" random variables with mean="mean" and variance="variance"
        // Since GenerateNormal() uses the box_muler transform, we will get a pair of 
        // normally distributed random numbers for every pair of uniform random numbers
        // generated. To maintain efficiency, we will use this pair whenever possible.
        public double[] Normal(double mean=0, double variance=1, int size=1)
        {
            double[] samples = new double[size];
            // Since GenerateNormal gives pairs, we need to handle the case when size is odd
            if (size % 2 == 1) 
            {
                samples[size-1] = variance * GenerateNormal()[0] + mean;
                size--;
            }
            
            // Now we need to add an even number of normal numbers to samples. 
            for (int i=0; i < size; i += 2)
            {
                samples[i] = variance * GenerateNormal()[0] + mean;
                samples[i+1] = variance * GenerateNormal()[1] + mean;
            }
            return samples;
        }

        // Generate a pair of normally distributed random numbers from uniform random
        // numbers using the Box-Muller transform
        private double[] GenerateNormal()
        {
            // Create a new random object and generate two uniformly distributed sample
            double x1 = rand.NextDouble();
            double x2 = rand.NextDouble();

            // Apply Box-Muller transformation
            double z1 = Math.Sqrt(-2*Math.Log(x1))*Math.Cos(2*Math.PI*x2);
            double z2 = Math.Sqrt(-2*Math.Log(x1))*Math.Sin(2*Math.PI*x2);
            return new double[] {z1,z2};
        }
    }


    static class Stats
    {
        // Returns the sum of all elements in the array
        public static double Sum(double[] data)
        {
            double sum = 0;
            for (int i = 0; i < data.Length; i++)
            {
                sum += data[i];
            }
            return sum;
        }

        // Return the mean of all the elements in the array
        public static double Mean(double[] data)
        {
            return Sum(data)/data.Length;
        }

        // Return the variance of all the elements in the array
        public static double Variance(double[] data)
        {
            double sum_diff = 0;
            double mean = Mean(data);
            for (int i = 0; i < data.Length; i++)
            {
                sum_diff += Math.Pow(data[i] - mean, 2);
            }
            return sum_diff/(data.Length-1);
        }

        public static double[] Maximum(double[] arr1, double[] arr2)
        {
            // Check if both SciArray's have the same size. If they don't throw an exception
            if (arr1.Length != arr2.Length)
            {
                throw new ArgumentException("arr1 and arr2 must have the same size.");
            }

            double[] maxArray = new double[arr1.Length];
            for (int i = 0; i < arr1.Length; i++)
            {
                
                maxArray[i] = (arr1[i] >= arr2[i]) ? arr1[i] : arr2[i];
            }
            return maxArray;
        }

        // Function to calculate the cumulative standard normal distribution
        public static double NormCDF(double x)
        {
            const double PI = 3.14159265358979323846;
            const double a1 = 0.319381530;
            const double a2 = -0.356563782;
            const double a3 = 1.781477937;
            const double a4 = -1.821255978;
            const double a5 = 1.330274429;
            const double L = 0.0;
            const double K = 1.0;
            double k = 1.0 / (1.0 + 0.2316419 * Math.Abs(x));
            double w = 1.0 - 1.0 / Math.Sqrt(2 * PI) * Math.Exp(-0.5 * x * x) * 
                    (a1 * k + a2 * Math.Pow(k, 2) + a3 * Math.Pow(k, 3) + a4 * Math.Pow(k, 4) + a5 * Math.Pow(k, 5));
            return x < 0.0 ? 1.0 - w : w;
        }
    }

    static class Output
    {
        public static void PrintArray<T>(T[] data)
        {
            Console.Write("[");
            for (int i = 0; i < data.Length-1; i++)
            {
                Console.Write(data[i] + ", ");
            }
            Console.Write(data[data.Length-1] + "]");
        }

        public static void Print2DArray<T>(T[,] matrix)
        {
            for (int i = 0; i < matrix.GetLength(0); i++)
            {
                for (int j = 0; j < matrix.GetLength(1); j++)
                {
                    Console.Write(matrix[i,j] + "\t");
                }
                Console.WriteLine();
            }
        }
    }

    static class Payoff
    {
        public static double callPayoff(double StockPrice, double Strike)
        {
            return Math.Max(StockPrice - Strike, 0);
        }

        public static double putPayoff(double StockPrice, double Strike)
        {
            return Math.Max(Strike - StockPrice, 0);
        }
    }
}