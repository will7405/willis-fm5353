// See https://aka.ms/new-console-template for more information

Thread t = new Thread(new ThreadStart(new Action() {
    Console.WriteLine("Hello World")
}));
t.start();
t.join();
