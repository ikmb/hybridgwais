/*
 *    Copyright (C) 2018-2023 by Lars Wienbrandt,
 *    Institute of Clinical Molecular Biology, Kiel University
 *    
 *    This file is part of HybridGWAIS.
 *
 *    HybridGWAIS is free software: you can redistribute it and/or modify
 *    it under the terms of the GNU General Public License as published by
 *    the Free Software Foundation, either version 3 of the License, or
 *    (at your option) any later version.
 *
 *    HybridGWAIS is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *    GNU General Public License for more details.
 *
 *    You should have received a copy of the GNU General Public License
 *    along with HybridGWAIS. If not, see <https://www.gnu.org/licenses/>.
 */

#ifndef STATUSFILE_H_
#define STATUSFILE_H_

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <stdio.h>

using namespace std;

class StatusFile {
public:
    static void updateFile(const string &statfile_, bool yaml) {
        instance.statfile.assign(statfile_);
        instance.infofile.assign(statfile_+".info");
        instance.warningfile.assign(statfile_+".warning");
        instance.errorfile.assign(statfile_+".error");
        instance.yaml = yaml;
        if (yaml) {
            instance.infofile.append(".yaml");
            instance.warningfile.append(".yaml");
            instance.errorfile.append(".yaml");
        }
        // clear files, if they already exist:
        // silent failure if the files do not exist or are not removable
        remove(instance.infofile.c_str());
        remove(instance.warningfile.c_str());
        remove(instance.errorfile.c_str());
    }

    static void updateStatus(float value, const string &info = "") {
        if (instance.statfile.empty())
            return;
        instance.statvalue = value;
        if (!info.empty())
            instance.statinfostr.assign(info);
        ofstream ofs(instance.statfile, ofstream::trunc);
        if (!ofs.fail())
            ofs << instance.statvalue << "," << instance.statinfostr << endl;
    }

    static void addInfoYAML(const string &key, const string &value) {
        stringstream s;
        if (instance.infoempty)
            s << "---\n"; // add YAML "header"
        s << "- " << key << ": " << value;
        if (!addToFile(instance.infofile, s.str()))
            addWarning("Could not write to info file.");
        else
            instance.infoempty = false;
    }

    static void addInfo(const string &info, bool yamlmessage = true) {
        cout << filterHTML(info) << endl;
        if ((instance.yaml && !yamlmessage ) || instance.statfile.empty())
            return;

        stringstream s;
        if (instance.infoempty && instance.yaml && yamlmessage)
            s << "---\n"; // add YAML "header"
        if (instance.yaml) {
            if (yamlmessage)
                addYAMLMessage(filterHTML(info), s);
        } else
            s << info << endl;
        if (!addToFile(instance.infofile, s.str()))
            addWarning("Could not write to info file.");
        else
            instance.infoempty = false;
    }

    static void addWarning(const string &warning) {
        cerr << "WARNING: " << filterHTML(warning) << endl;
        if (instance.statfile.empty())
            return;

        stringstream s;
        if (instance.warningempty) {
            if (instance.yaml)
                s << "---\n"; // add YAML "header"
            else
                s << "<h3>WARNING:</h3>\n";
        }
        if (instance.yaml)
            addYAMLMessage(filterHTML(warning), s);
        else
            s << "<ul><li>" << warning << "</li></ul>" << endl;
        if (!addToFile(instance.warningfile, s.str()))
            cerr << "WARNING: Could not write to warning file." << endl;
        else
            instance.warningempty = false;

        if (!instance.yaml) { // warning messages will not appear in YAML info file
            stringstream si;
            si << "<b>WARNING:</b> " << warning << "<br>" << endl;
            if (!addToFile(instance.infofile, si.str()))
                cerr << "WARNING: Could not write to info file." << endl;
            else
                instance.infoempty = false;
        }
    }

    static void addError(const string &error) {
        cerr << "ERROR: " << filterHTML(error) << endl;
        if (instance.statfile.empty())
            return;

        stringstream s;
        if (instance.errorempty) {
            if (instance.yaml)
                s << "---\n"; // add YAML "header"
            else
                s << "<h3>ERROR:</h3>\n";
        }
        if (instance.yaml)
            addYAMLMessage(filterHTML(error), s);
        else
            s << "<ul><li>" << error << "</li></ul>" << endl;
        if (!addToFile(instance.errorfile, s.str()))
            addWarning("Could not write to error file.");
        else
            instance.errorempty = false;

        if (!instance.yaml) { // error messages will not appear in YAML info file
            stringstream si;
            si << "<b>ERROR:</b> " << error << "<br>" << endl;
            if (!addToFile(instance.infofile, si.str()))
                addWarning("Could not write to info file.");
            else
                instance.infoempty = false;
        }

        // Automatically set status to "Failed" with current progress value + string
        updateStatus(instance.statvalue, "Failed: " + instance.statinfostr);
    }


private:
    StatusFile(){}

    // returns true on success, false on failure
    static bool addToFile(const string &filename, const string &s) {
        ofstream ofs(filename, ofstream::app);
        if (!ofs.fail()) {
            ofs << s << endl;
            return true;
        } else
            return false;
    }

    // simple filter for HTML Tags <>
    static string filterHTML(const string &s_) {
        string s(s_);
        size_t end = s.rfind('>'); // find a closing bracket first
        while (end != string::npos) {
            size_t start = s.rfind('<', end); // find corresponding opening bracket
            if (start == string::npos) // didn't find an opening bracket
                break;
            s.erase(start, end-start+1);
            end = s.rfind('>', start);
        }
        return s;
    }

    static void addYAMLMessage(const string &message, stringstream &s) {
        const string ind("    "); // indentation
        s << "- message:\n"; // message key
        // find newline characters and insert indentation for each one
        size_t start = 0, pos = 0;
        while(pos != string::npos) {
            pos = message.find('\n', start);
            size_t spos = message.find_first_not_of(" \t", start); // find position after leading whitespace
            if (spos == string::npos)
                break;
            string m = message.substr(spos, pos - spos);
            if (!m.empty()) // only print non-empty messages
                s << ind << m << "\n";
            start = pos+1; // character after newline
        }
    }

    string statfile;
    string statinfostr;
    float statvalue = 0.0;

    string infofile;
    string warningfile;
    string errorfile;

    bool infoempty = true;
    bool warningempty = true;
    bool errorempty = true;

    bool yaml = false;

    static StatusFile instance;
};

#endif /* STATUSFILE_H_ */
