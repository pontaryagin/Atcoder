using System.Text;
using System.Net;
using System.Text.RegularExpressions;
using static System.Console;
using System.IO;

namespace getTestcaseForAtcoder
{
	public class Program
	{
		public static void Main(string[] args)
		{
			WriteLine("Input Atcoder Problem URL : ");
			string testcases = getTestCase(getContext(ReadLine().Trim(new char[] { '\r', '\n', ' ' })));

			using (StreamWriter sw = new StreamWriter("ProgramTests.cs"))
			{
				sw.Write(testcases);
			}

			WriteLine("Output to \"ProgramTests.cs\"");
		}

		static string getTestCase(string html)
		{
			string source1 =
@"using Microsoft.VisualStudio.TestTools.UnitTesting;
using System;
using System.IO;

namespace Atcoder.Tests
{
    [TestClass()]
    public class ProgramTests
    {
        (TESTCASES)
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
}";

			string source2 = @"
        [TestMethod()]
        [Timeout(2000)]
        public void Test(TESTNUM)()
        {
            string input  = ""(TESTCASE1)"";
            string output = ""(TESTCASE2)"";

            AssertIO(input, output);
        }";

			string anchor = "(<pre class=\"prettyprint linenums\">|<pre>)(?<testcase>.*?)</pre>";
			Regex re = new Regex(anchor, RegexOptions.IgnoreCase | RegexOptions.Singleline);

			StringBuilder sb = new StringBuilder();
			int i = 1;
			for (Match m = re.Match(html); m.Success; m = m.NextMatch())
			{
				if (m.Groups["testcase"].Value.Contains("<var>")) continue;

				string testCase1 = m.Groups["testcase"].Value; m = m.NextMatch();
				string testCase2 = m.Groups["testcase"].Value;

				testCase1 = testCase1.TrimStart('\r').TrimStart('\n').Replace("\r", "").Replace("\n", "\\n");
				testCase2 = testCase2.TrimStart('\r').TrimStart('\n').TrimEnd('\n').TrimEnd('\r').Replace("\r", "").Replace("\n", "\\n");

				sb.AppendLine(source2.Replace("(TESTNUM)", (i++).ToString()).Replace("(TESTCASE1)", testCase1).Replace("(TESTCASE2)", testCase2));
			}

			return source1.Replace("(TESTCASES)", sb.ToString());
		}

		static string getContext(string url)
		{
			WebClient wc = new WebClient();
			byte[] data = wc.DownloadData(url);
			Encoding enc = Encoding.GetEncoding("utf-8");
			string ret = enc.GetString(data);
			wc.Dispose();
			return ret;
		}
	}
}