%> @file classDocumentationExample.m
%> @brief File used to show an example of class description
% ======================================================================
%> @brief Here we have a brief description of the class.
%
%> And here we can put some more detailed informations about the class.
% ======================================================================
classdef  (InferiorClasses = {?class1,?class2}) classDocumentationExample

  properties (Access = protected)
    %> Description of a protected property
    protectedProperty
  end
  properties (Access = public)
    %> Description of a public property
    publicProperty
  end
  properties (Access = private)
    %> Description of a private property
    pivateProperty
  end
  properties (Constant = true)
    %> Description of a constant property
    constantProperty = {'1', '2', 'trois'};
  end
  properties
    %> Description of the first property of the class
    first_property = []
    %> Description of the second property of the class
    second_property = []
    %> Description of the third property of the class
    third_property = [1  2];
  end
  events
    %> Description of first event
    FirstEvent
    %> Description of second event
    SecondEvent
  end
  %> Description of the enumeration.
  enumeration
    %> Description of the first item
    one (1)
    %> Description of the second item
    two (2)
    %> Description of the third item
    three
  end
  methods
    % ======================================================================
    %> @brief Class constructor
    %>
    %> More detailed description of what the constructor does.
    %>
    %> @param param1 Description of first parameter
    %> @param anotherParam Description of the second parametere
    %>
    %> @return instance of the classDocumentationExample class.
    % ======================================================================
    function obj = classDocumentationExample(param1, anotherParam)
    end

    % ======================================================================
    %> @brief Brief description of the exampleMethod1 method
    %>
    %> @param obj instance of the classDocumentationExample class.
    % ======================================================================
    function exampleMethod1(obj)
    end

    % ======================================================================
    %> @brief Brief description of the exampleMethod2 method
    %>
    %> @param obj instance of the classDocumentationExample class.
    %> @retval ret return value of this method
    % ======================================================================
    function ret = exampleMethod2(obj)
    end
  end
  methods (Static=true)
    % ======================================================================
    %> @brief Brief description of the exampleStaticMethod method
    %>
    %> This method is static and public, with an inused (~) argument
    %> @param param1 Description of the parameter
    %> @param param2 Description of the parameter
    %> @retval out return value of this method
    % ======================================================================
    function out = exampleStaticMethod(param1, ~, param2)
    end
  end
  methods (Static, Access=private)
    % ======================================================================
    %> @brief Brief description of the exampleStaticPrivateMethod method
    %>
    %> This method is static and private
    %> @param param1 Description of the parameter
    %> @param param2 Description of the parameter
    %> @retval out return value of this method
    % ======================================================================
    function out = exampleStaticPrivateMethod(param1, param2)
    end
  end
  methods (Access=protected, Static)
    % ======================================================================
    %> @brief Brief description of the exampleStaticProtectedMethod method
    %>
    %> This method is static and protected
    %> @param param1 Description of the parameter
    %> @retval out return value of this method
    % ======================================================================
    function out = exampleStaticProtectedMethod(param1)
    end
  end
  methods (Access=public, Static = true)
    % ======================================================================
    %> @brief Brief description of the exampleStaticPublicMethod method
    %>
    %> This method is static and public
    %> @param param1 Description of the parameter
    %> @retval out return value of this method
    % ======================================================================
    function out = exampleStaticPublicMethod(param1)
    end
  end
  methods (Access=private)
    % ======================================================================
    %> @brief Brief description of the examplePrivateMethod method
    %>
    %> This method is private
    %> @param param1 Description of the parameter
    %> @retval out return value of this method
    % ======================================================================
    function out = examplePrivateMethod(param1)
    end
  end
  methods (Access=protected)
    % ======================================================================
    %> @brief Brief description of the exampleProtectedMethod method
    %>
    %> This method is protected
    %> @param param1 Description of the parameter
    %> @retval out return value of this method
    % ======================================================================
    function out = exampleProtectedMethod(param1)
    end
  end
  methods (Access=public, Static=false)
    % ======================================================================
    %> @brief Brief description of the examplePublicMethod2 method
    %>
    %> This method is public and not static
    %> @param param1 Description of the parameter
    %> @retval out return value of this method
    % ======================================================================
    function out = examplePublicMethod2(param1)
    end
  end
  methods (Access=public, ~Static)
    % ======================================================================
    %> @brief Brief description of the exampleNonStaticPublicMethod3 method
    %>
    %> This method is public and not static
    %> @param param1 Description of the parameter
    %> @retval out return value of this method
    % ======================================================================
    function out = exampleNonStaticPublicMethod3(param1)
    end
  end
  methods (Abstract = true)
    % ======================================================================
    %> @brief Brief description of the exampleAbstractMethod method
    %>
    %> This method is abstract : only the signature of this function is
    %> declared.
    %> @param param1 Description of the first parameter
    %> @param param2 Description of the second parameter
    %> @retval out return value of this method
    % ======================================================================
    out = exampleAbstractMethod(param1, param2);
end
