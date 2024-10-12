using Microsoft.EntityFrameworkCore;
using Microsoft.EntityFrameworkCore.Design;
using System;
using System.ComponentModel.DataAnnotations.Schema;

namespace MyApp
{
    internal class Program
    {
        static void Main(string[] args)
        {
            
        }
    }

    public class Student
    {
        public int Id { get; set; }
        [Column("studentname")]
        public string Name { get; set; }
    }
}