/**
@file Exceptions.h
@brief Declaration of the various exception classes.

THAMES tries to implement exception handling consistently throughout the code,
although more could be done to check for out-of-bounds errors on arrays and
other container classes.

Right now, many of these classes are quite similar and could probably be structured
better with a generic base class, but the thought is that each class will take
on more distinct details later.

@todo Make a virtual base class for common things like descriptions, etc.

*/

#ifndef EXCEPTIONSH
#define EXCEPTIONSH

#include <iostream>
#include <string>

using namespace std;

/**
@class Declare the EOBException class

The `EOBException` class handles exceptions caused by attempting to access
an out of bounds element of an array.
*/
class EOBException {

private:

    string arrayname_;      /**< Name of the array accessed */
    string classname_;      /**< Name of the class that accessed the array */
    string functionname_;   /**< Name of the method that accessed the array */
    int sizelimit_;         /**< Number of elements contained in the array */
    unsigned int indx_;     /**< Out-of-bounds element number that was queried */

public:

/**
@brief Default constructor initializes class members to default (blank) values.

*/
EOBException ()
{
    classname_ = "";
    functionname_ = "";
    arrayname_ = "";
    sizelimit_ = 0;
    indx_ = 0;
}

/**
@brief Overloaded constructor that is typically invoked by THAMES.

@param cname is the class name where the exception was thrown
@param fname is the method name where the exception was thrown
@param arname is the name of the array that was accessed
@param sl is the total number of array elements in that array
@param id is the element number (out of bounds) that was queried erroneously
*/
EOBException (const string &cname,
              const string &fname,
              const string &arname,
              const int sl,
              const unsigned int id)
{
    classname_ = cname;
    functionname_ = fname;
    arrayname_ = arname;
    sizelimit_ = sl;
    indx_ = id;
}

/**
@brief Get the class name responsible for throwing the EOB exception.

@return the class name
*/
string &getClassname () const
{
    return (string &)classname_;
}

/**
@brief Get the function name responsible for throwing the EOB exception.

@return the function name
*/
string &getFunctionname () const
{
    return (string &)functionname_;
}

/**
@brief Get the array name that was queried when the EOB exception was thrown.

@return the array name
*/
string &getArrayname () const
{
    return (string &)arrayname_;
}

/**
@brief Get the number of elements of the array queried when the EOB exception was thrown.

@return the array name
*/
int getSizelimit() const { return sizelimit_; }

/**
@brief Get the queried (out-of-bounds) index number.

@return the queried index number
*/
unsigned int getIndx () const
{
    return indx_;
}

/**
@brief Provide formatted output of the exception details.

*/
void printException ()
{
        cout << endl << "EOB Exception Thrown:" << endl;
        cout << "    Details: " << endl;
        cout << "        Offending Function " << classname_ << "::" << functionname_ << endl;
        cerr << endl << "EOB Exception Thrown:" << endl;
        cerr << "    Details: " << endl;
        cerr << "        Offending Function " << classname_ << "::" << functionname_ << endl;
        if (indx_ == 0) {
            cout << "        Array: " << arrayname_ << endl;
            cerr << "        Array: " << arrayname_ << endl;
        } else {
            cout << "        Array: " << arrayname_ << " contains " << sizelimit_;
            cout << " elements, but tried to access element " << indx_ << endl;
            cerr << "        Array: " << arrayname_ << " contains " << sizelimit_;
            cerr << " elements, but tried to access element " << indx_ << endl;
        }
        return;
}

};  // End of the EOBException class


/**
@class Declare the FileException class

The `FileException` class handles exceptions caused by attempting to open, close,
write to, or read from a file.

*/
class FileException {

private:

    string filename_;       /**< Name of the offending file */
    string extype_;         /**< Name of the exception description */
    string classname_;      /**< Name of the class that threw the exception */
    string functionname_;   /**< Number of function that threw the exception */

public:

/**
@brief Default constructor initializes class members to default (blank) values.

*/
FileException ()
{
    classname_ = "";
    functionname_ = "";
    filename_ = "";
    extype_ = "";
}

/**
@brief Overloaded constructor that is typically invoked by THAMES.

@param cname is the class name where the exception was thrown
@param fname is the method name where the exception was thrown
@param filename is the name of the offending file
@param extype is the description of the exception type
*/
FileException (const string &cname,
               const string &fname,
               const string &filename,
               const string &extype)
{
    classname_ = cname;
    functionname_ = fname;
    filename_ = filename;
    extype_ = extype;
}

/**
@brief Get the class name responsible for throwing the file exception.

@return the class name
*/
string &getClassname () const
{
    return (string &)classname_;
}

/**
@brief Get the function name responsible for throwing the file exception.

@return the function name
*/
string &getFunctionname () const
{
    return (string &)functionname_;
}

/**
@brief Get the file name that was queried when the file exception was thrown.

@return the file name
*/
string &getFilename () const
{
    return (string &)filename_;
}

/**
@brief Get the file exception type description.

@return the exception type description
*/
string &getExtype () const
{
    return (string &)extype_;
}

/**
@brief Provide formatted output of the exception details.

*/
void printException()
{
        cout << endl << "File Exception Thrown:" << endl;
        cout << "    Details: " << endl;
        cout << "        Offending Function " << classname_ << "::" << functionname_ << endl;
        cout << "        File: " << filename_ << ", Problem:" << extype_ << endl;
        cerr << endl << "GEM Exception Thrown:" << endl;
        cerr << "    Details: " << endl;
        cerr << "        Offending Function " << classname_ << "::" << functionname_ << endl;
        cerr << "        File: " << filename_ << ", Problem: " << extype_ << endl;
        return;
}

};      // End of the FileException class

/**
@class Declare the FloatException class

The `FloatException` class handles exceptions caused by floating point operations,
especially divide-by-zero exceptions.

*/
class FloatException {

private:

    string description_;    /**< Description of the floating point exception */
    string classname_;      /**< Name of the class that threw the exception */
    string functionname_;   /**< Number of function that threw the exception */

public:

/**
@brief Default constructor initializes class members to default (blank) values.

*/
FloatException ()
{
    classname_ = "";
    functionname_ = "";
    description_ = "";
}

/**
@brief Overloaded constructor that is typically invoked by THAMES.

@param cname is the class name where the exception was thrown
@param fname is the method name where the exception was thrown
@param strd is the description of the exception
*/
FloatException (const string &cname,
                const string &fname,
                const string &strd)
{
    classname_ = cname;
    functionname_ = fname;
    description_ = strd;
}

/**
@brief Get the class name responsible for throwing the floating point exception.

@return the class name
*/
string &getClassname () const
{
    return (string &)classname_;
}

/**
@brief Get the function name responsible for throwing the floating point exception.

@return the function name
*/
string &getFunctionname () const
{
    return (string &)functionname_;
}

/**
@brief Get the description of the floating point exception.

@return the file name
*/
string &getDescription () const
{
    return (string &)description_;
}

/**
@brief Provide formatted output of the exception details.

*/
void printException ()
{
        cout << endl << "Floating Point Exception Thrown:" << endl;
        cout << "    Details: " << endl;
        cout << "        Offending Function " << classname_ << "::" << functionname_ << endl;
        cout << "        Description: " << description_ << endl;
        cerr << endl << "Floating Point Exception Thrown:" << endl;
        cerr << "    Details: " << endl;
        cerr << "        Offending Function " << classname_ << "::" << functionname_ << endl;
        cerr << "        Description: " << description_ << endl;
        return;
}

};      // End of the FloatException class


/**
@class Declare the HandleException class

The `HandleException` class handles exceptions caused by errors in dealing
with data handles.

*/
class HandleException {

private:

    string description_;    /**< Description of the handle exception */
    string classname_;      /**< Name of the class that threw the exception */
    string functionname_;   /**< Number of function that threw the exception */
    string handle_;         /**< Description of the handle causing the exception */

public:

/**
@brief Default constructor initializes class members to default (blank) values.

*/
HandleException ()
{
    classname_ = "";
    functionname_ = "";
    handle_ = "";
    description_ = "";
}

/**
@brief Overloaded constructor that is typically invoked by THAMES.

@param cname is the class name where the exception was thrown
@param fname is the method name where the exception was thrown
@param handle is the handle that caused the exception
@param strd is the description of the exception
*/
HandleException (const string &cname,
                 const string &fname,
                 const string &handle,
                 const string &strd)
{
    classname_ = cname;
    functionname_ = fname;
    handle_ = handle;
    description_ = strd;
}

/**
@brief Get the class name responsible for throwing the handle exception.

@return the class name
*/
string &getClassname () const
{
    return (string &)classname_;
}

/**
@brief Get the function name responsible for throwing the handle exception.

@return the function name
*/
string &getFunctionname () const
{
    return (string &)functionname_;
}

/**
@brief Get the handle causing the exception.

@return the handle
*/
string &getHandle () const
{
    return (string &)handle_;
}

/**
@brief Get the description of the handle exception.

@return the file name
*/
string &getDescription () const
{
    return (string &)description_;
}

/**
@brief Provide formatted output of the exception details.

*/
void printException ()
{
    cout << endl << "Handle Exception Thrown:" << endl;
    cout << "    Details: " << endl;
    cout << "        Offending Function " << classname_ << "::" << functionname_ << endl;
    cout << "        Description: " << description_ << endl;
    cout << "             Handle: " << handle_ << endl;
    cerr << endl << "Floating Point Exception Thrown:" << endl;
    cerr << "    Details: " << endl;
    cerr << "        Offending Function " << classname_ << "::" << functionname_ << endl;
    cerr << "        Description: " << description_ << endl;
    cerr << "             Handle: " << handle_ << endl;
    return;
}

};      // End of the HandleException class


/**
@class Declare the GEMException class

The `GEMException` class handles exceptions caused by errors originating
in the GEM3K library.

*/
class GEMException {

private:

    string description_;    /**< Description of the GEM exception */
    string classname_;      /**< Name of the class that threw the exception */
    string functionname_;   /**< Number of function that threw the exception */

public:

/**
@brief Default constructor initializes class members to default (blank) values.

*/
GEMException ()
{
    classname_ = "";
    functionname_ = "";
    description_ = "";
}

/**
@brief Overloaded constructor that is typically invoked by THAMES.

@param cname is the class name where the exception was thrown
@param fname is the method name where the exception was thrown
@param strd is the description of the exception
*/
GEMException (const string &cname,
              const string &fname,
              const string &strd)
{
    classname_ = cname;
    functionname_ = fname;
    description_ = strd;
}

/**
@brief Get the class name responsible for throwing the handle exception.

@return the class name
*/
string &getClassname () const
{
    return (string &)classname_;
}

/**
@brief Get the function name responsible for throwing the handle exception.

@return the function name
*/
string &getFunctionname () const
{
    return (string &)functionname_;
}

/**
@brief Get the description of the handle exception.

@return the file name
*/
string &getDescription () const
{
    return (string &)description_;
}

/**
@brief Provide formatted output of the exception details.

*/
void printException ()
{
    cout << endl << "GEM Exception Thrown:" << endl;
    cout << "    Details: " << endl;
    cout << "        Offending Function " << classname_ << "::" << functionname_ << endl;
    cout << "        " << description_ << endl;
    cerr << endl << "GEM Exception Thrown:" << endl;
    cerr << "    Details: " << endl;
    cerr << "        Offending Function " << classname_ << "::" << functionname_ << endl;
    cerr << "        " << description_ << endl;
    return;
}

};      // End of GEMException class


/**
@class Declare the DataException class

The `DataException` class handles exceptions related to miscellaneous
data errors.

*/
class DataException {

private:

    string description_;    /**< Description of the GEM exception */
    string classname_;      /**< Name of the class that threw the exception */
    string functionname_;   /**< Number of function that threw the exception */

public:

/**
@brief Default constructor initializes class members to default (blank) values.

*/
DataException ()
{
    classname_ = "";
    functionname_ = "";
    description_ = "";
}

/**
@brief Overloaded constructor that is typically invoked by THAMES.

@param cname is the class name where the exception was thrown
@param fname is the method name where the exception was thrown
@param strd is the description of the exception
*/
DataException (const string &cname,
               const string &fname,
               const string &strd)
{
    classname_ = cname;
    functionname_ = fname;
    description_ = strd;
}

/**
@brief Get the class name responsible for throwing the handle exception.

@return the class name
*/
string &getClassname () const
{
    return (string &)classname_;
}

/**
@brief Get the function name responsible for throwing the handle exception.

@return the function name
*/
string &getFunctionname () const
{
    return (string &)functionname_;
}

/**
@brief Get the description of the handle exception.

@return the file name
*/
string &getDescription () const
{
    return (string &)description_;
}

/**
@brief Provide formatted output of the exception details.

*/
void printException ()
{
    cout << endl << "Data Exception Thrown:" << endl;
    cout << "    Details: " << endl;
    cout << "        Offending Function " << classname_ << "::" << functionname_ << endl;
    cout << "        Problem:" << description_ << endl;
    cerr << endl << "Data Exception Thrown:" << endl;
    cerr << "    Details: " << endl;
    cerr << "        Offending Function " << classname_ << "::" << functionname_ << endl;
    cerr << "        Problem: " << description_ << endl;
    return;
}

};      // End of DataException class

#endif
