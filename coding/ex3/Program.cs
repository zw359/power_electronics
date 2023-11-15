// See https://aka.ms/new-console-template for more information
//Console.WriteLine("Hello, World!");
using System;
using System.Xml;
namespace ex3
{
	static class Program
	{
		static void Main()
		{
			XmlDocument doc = new XmlDocument();
			doc.Load("test2.xml");
			
			XmlNamespaceManager nsmgr = new XmlNamespaceManager(doc.NameTable);
			nsmgr.AddNamespace("df","urn:schemas-microsoft-com:office:spreadsheet");
			nsmgr.AddNamespace("o","urn:schemas-microsoft-com:office:office");
			nsmgr.AddNamespace("x","urn:schemas-microsoft-com:office:excel");
			nsmgr.AddNamespace("dt","uuid:C2F41010-65B3-11d1-A29F-00AA00C14882");
			nsmgr.AddNamespace("ss","urn:schemas-microsoft-com:office:spreadsheet");
			nsmgr.AddNamespace("html","http://www.w3.org/TR/REC-html40");
			//XmlNode node = doc.SelectSingleNode("/df:Workbook/df:Worksheet", nsmgr);
			XmlNode node = doc.SelectSingleNode("/df:Workbook/df:Worksheet[@ss:Name='Input Products']", nsmgr);
			//XmlNode node = doc.SelectSingleNode("/df:Workbook/df:Worksheet[last()]", nsmgr);
			//XmlNode node = doc.SelectSingleNode("/df:Workbook/df:Worksheet[last()]/df:Table/df:Row[last()]", nsmgr);
			//XmlNode node = doc.SelectSingleNode("/df:Workbook/df:Worksheet[@ss:Name='Input Products']/df:Table/df:Row[8]/df:Cell/df:Data[1]", nsmgr);
			
			//Console.WriteLine(node.Attributes[0].InnerText);
			//Console.WriteLine(node.Attributes[0].Name);
			//Console.WriteLine(node.Name);
			//Console.WriteLine(node.InnerText);
			Console.WriteLine(node.GetAttribute("ss:Name", nsmgr));
			//Console.WriteLine(doc.DocumentElement.ChildNodes[0].ChildNodes[0].Name);
		
			
		}
	}
}
