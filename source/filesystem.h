#ifndef FILESYSTEM_H
#define	FILESYSTEM_H

namespace FS
{
	static string root = "D:/ResultsRaw/sim/";
	static string result = root + "result/";
	static string image = root + "image/";
	static string save = root + "save/";

	static void setFS(Setup &s)
	{
		root = s.opt.outputPath + "/" + s.run.name + "/";

		result = root + "result/";
		image = root + "image/";
		save = root + "save/";
	}

	static void createFS()
	{
		CreateDirectory (root.c_str(), NULL);
		CreateDirectory (result.c_str(), NULL);
		CreateDirectory (image.c_str(), NULL);
		CreateDirectory (save.c_str(), NULL);
	}

	static string getFN(string dir, string name, string ext, int num = -1)
	{
		ostringstream ss;
		if(num < 0)
			ss << dir << name << "." << ext << flush;
		else
			ss << dir << name << "_" << setfill('0') << setw(12) << num << "." << ext << flush;
		return ss.str();
	}

	static bool dirExists(const string &dirName_in)
	{
		DWORD ftyp = GetFileAttributesA(dirName_in.c_str());
		if (ftyp == INVALID_FILE_ATTRIBUTES)
			return false;  //something is wrong with your path!

		if (ftyp & FILE_ATTRIBUTE_DIRECTORY)
			return true;   // this is a directory!

		return false;    // this is not a directory!
	}

	static bool dirCreate(const string &dirName_in)
	{
		if(CreateDirectoryA(dirName_in.c_str(), NULL))
		{
			string msg = "Folder created: " + dirName_in;
			MessageBoxA(NULL, msg.c_str(), "", NULL);
			return true;
		}
		else
		{
			string msg = "Failed to create folder: " + dirName_in;
			MessageBoxA(NULL, msg.c_str(), "", NULL);
			return false;
		}
	}
}



#endif