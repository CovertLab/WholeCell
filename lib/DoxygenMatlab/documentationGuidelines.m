%> @file  documentationGuidelines.m
%> @brief Documentation guidelines to extract documentation using Doxygen
%======================================================================
%> @mainpage Documentation guidelines
%>
%> @section intro Introduction
%>
%> The @b doxygen software (http://www.doxygen.org) allows you to extract code comments
%> from your code to automatically generate documentation.
%>
%> Doxygen builds automatically a documentation based on C++ components (classes,
%> functions, structs, variables,...) using syntaxic analysis. Doxygen extends this
%> documentation by extracting special comments.
%>
%> @section matlabComments Extract code comments from Matlab files (.m files)
%>
%> Matlab code is not supported nativly by Doxygen : a perl filter allowing the conversion
%> from .m to C++ has been developped, the @c m2cpp.pl script.
%>
%> This filter extracts only :
%> - lines beginning with <b>%></b> :
%>   these lines are converted into C++ comments, ie beginning by ///
%> - lines beginning with the @c function keyword :
%>   these lines are converted into C++ functions
%> - lines beginning with the @c classdef keyword :
%>   these lines are converted into C++ classes
%> - lines beginning with the @c properties keyword :
%>   these lines are converted into C++ properties
%> - lines beginning with the @c enumeration keyword :
%>   these lines are converted into C++ enum
%> - lines beginning with the @c events keyword :
%>   these lines are converted into C++ enum Events
%>   @note for enumeration definition : this script only supports the following declaration :
%>   @code
%>   enumeration
%>     first_enum
%>     second_enum
%>   end
%>   @endcode
%>   This one is not yet supported :
%>   @code
%>   enumeration
%>     first_enum, second_enum
%>   end
%>   @endcode
%>   The file classDocumentationExample.m provides an example of enumeration definition
%>   with comments extracted by Doxygen.
%>
%> .
%>
%> See @ref installation if you want more details about how to make this script work.
%>
%> @attention Each line belonging to the doxygen documentation must begin with <b>%></b> .
%>
%> @par Example
%>
%>@code
%>% Matlab comment ignored by doxygen
%>%> comment analyzed by doxygen
%>@endcode
%> @attention Doxygen keyword have to begin by @b @, for example @@b to bold the text (the use of \ instead of @ is not supported)
%>
%>@section funcDecr Function description
%> The keyword @b @@param and @b @@retval will be used to describe the input and
%> output parameters of a function.
%>
%> For function description, the description should follow the following presentation :
%>
%>@verbatim
%>% ======================================================================
%>%> @brief Brief description of the function
%>%>
%>%> More detailed description.
%>%>
%>%> @param arg1 First argument
%>%> @param arg2 Second argument
%>%>
%>%> @retval out1 return value for the first output variable
%>%> @retval out2 return value for the second output variable
%>% ======================================================================
%>[out1, out2] = function( arg1, arg2)
%>  out1 = arg2;
%>  out2 = arg1;
%>end
%>@endverbatim
%>
%>
%> @section classDecr Class description
%>
%> For class description, the following description can be used :
%>
%>@verbatim
%>% ======================================================================
%>%> @brief Brief description of the class
%>%>
%>%> More detailed description of the class
%>%>
%>% ======================================================================
%>@endverbatim
%>
%>
%> @page installation Installation details
%>
%> This package contains two files :
%> - a perl script (@c m2cpp.pl) which is a filter that converts  .m files into .cpp files that Doxygen can understand
%> - the @c Doxyfile file (configuration file used by Doxygen) that contains parameters needed by Doxygen to extract the documentation from .m files.
%>
%> Installation :
%> - You need to have the @b Doxygen software installed (version 1.5.9 or newer required (tested with version 1.7.1))
%> - You need to have @b perl installed (perl is shipped with Matlab, located usually  in @c $matlabroot\\sys\\perl\\win32\\bin)
%> - extract the @c Doxyfile file from the @c doxygen.zip package and replace the default Doxyfile provided by Doxygen
%> - extract the @c m2cpp.pl into a directory (for example @c C:\\DoxygenMatlab)
%> - edit the Doxyfile file (or use the DoxyWizard tool provided by Doxygen) to modify a few settings :
%>   - EXTENSION_MAPPING=.m=C++
%>   - FILTER_PATTERN=*m=C:\\DoxygenMatlab\\m2cpp.pl
%>   - PERL_PATH=\<path to your perl version\>
%>   - INPUT=\<directory where your documented code is located\>
%>   - OUTPUT_DIRECTORY=\<directory where you want to generate your documentation\>
%>   - STRIP_FORM_PATH=\<directory where your documented code is located\>
%>
%>   @note For Windows users, in certain circumstances, the association
%>   between .pl files and the perl executable is not well configured,
%>   leading to "Argument must contain filename -1 at C:\DoxygenMatlab\m2cpp.pl line 4"
%>   when running doxygen. To work around this issue, you should execute the following lines
%>   in a Windows command prompt ("cmd") :
%>   @code
%>   assoc .pl=PerlScript
%>   ftype PerlScript=C:\Program Files\MATLAB\R2010b\sys\perl\win32\bin\perl.exe %1 %*
%>   @endcode
%>   (don't forget to replace the path to the perl.exe file with yours in the line above)
%>
%>
%>

%======================================================================
%> @brief Brief description of the function
%>
%> More detailed description.
%>
%> @param arg1 First argument
%> @param arg2 Second argument
%>
%> @retval out1 return value for the first output variable
%> @retval out2 return value for the second output variable
%======================================================================
function [out1, out2] = test(arg1, arg2)
end
