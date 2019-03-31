using Microsoft.VisualStudio.TestTools.UnitTesting;
using System;
using System.IO;

namespace Atcoder.Tests
{
    [TestClass()]
    public class ProgramTests
    {
        
        [TestMethod()]
        [Timeout(2000)]
        public void Test1()
        {
            string input  = "A\n";
            string output = "T";

            AssertIO(input, output);
        }

        [TestMethod()]
        [Timeout(2000)]
        public void Test2()
        {
            string input  = "G\n";
            string output = "C";

            AssertIO(input, output);
        }

        [TestMethod()]
        [Timeout(2000)]
        public void Test3()
        {
            string input  = "A\n";
            string output = "T";

            AssertIO(input, output);
        }

        [TestMethod()]
        [Timeout(2000)]
        public void Test4()
        {
            string input  = "G\n";
            string output = "C";

            AssertIO(input, output);
        }

        private void AssertIO(string input, string output)
        {
            StringReader reader = new StringReader(input);
            Console.SetIn(reader);

            StringWriter writer = new StringWriter();
            Console.SetOut(writer);

            Program.Main(new string[0]);

            Assert.AreEqual(output + Environment.NewLine, writer.ToString());
        }
    }
}