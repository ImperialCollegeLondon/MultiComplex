classdef multicomplex % Of the form zn:[a1,a2,...,an]
    %{
    Copyright 2018, 2019 Jose Maria Varas Casado and Robert Hewson,
    Imperial College London.
                            
    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.
    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.
    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.
    %}
    % edited by Chung Chung To in 2021
    %% Properties %%
    
    properties
        zn         % The object zn is a matix of the form zn:[a1,a2,...,an]
    end
    
     % Flag for check of convergence and check for zero divisors -
     % 'F' : False, 'T' = True
     % Recomended to leave False
     properties (Constant = true)
       check_convergence = 'F';
       check_Z0 = 'F';
    end
    
    %% Initialization %%
    
    methods
        function self=multicomplex(a)         % Returns error if there are more than one input
            if nargin ~= 1
                error('Requires 1 array input')
            end
            if mod(log2(length(a)),1) ~= 0    % Multicomplex numbers need 2^n coefficients, for C_n
                error('Requires 2^n input size')
            end
            
            self.zn=a;                         % Initialization of ClassDef
            
        end
    end
    
   
    methods
  %% Basic Operator Overloading %%
  
        function mrep = matrep(self)                              % Matrix Representation
            
            % The matrix representation can only go up to a certain C_n (in
            % this case C_10). You could change this easily by adding more if
            % statements inside the while loop. However, MATLAB has a limit in
            % the memory used to store a matrix, so for n>13, the number of
            % elements inside the matrix would more than 2^14 x 2^14, which is more
            % than the maximum allowed. This imposes an upper limit on the
            % order of the multicomplex number that MATLAB can perform
            % operations with.
            
            sz=length(self.zn);
            
            n=log2(sz);
            
            if length(self.zn) == 1
                
                mrep=[self.zn];
                
            end
            
            i=2;
            
            while i<=sz
                
                m1=self.zn(i-1);
                m2=self.zn(i);
                
                C1=[m1,-m2;m2,m1];
                if rem(i,2^2) == 0
                    C2 = [storeC1,-C1;C1,storeC1];
                if rem(i,2^3) == 0
                    C3 = [storeC2,-C2;C2,storeC2];
                if rem(i,2^4) == 0
                    C4 = [storeC3,-C3;C3,storeC3];
                if rem(i,2^5) == 0
                    C5 = [storeC4,-C4;C4,storeC4];
                if rem(i,2^6) == 0
                    C6 = [storeC5,-C5;C5,storeC5];
                if rem(i,2^7) == 0
                    C7 = [storeC6,-C6;C6,storeC6];
                if rem(i,2^8) == 0
                    C8 = [storeC7,-C7;C7,storeC7];
                if rem(i,2^9) == 0
                    C9 = [storeC8,-C8;C8,storeC8];
                if rem(i,2^10) == 0
                    C10 = [storeC9,-C9;C9,storeC9];
                if rem(i,2^11) == 0
                    C11 = [storeC10,-C10;C10,storeC10];
                if rem(i,2^12) == 0
                    C12 = [storeC11,-C11;C11,storeC11];
                if rem(i,2^13) == 0
                    C13 = [storeC12,-C12;C12,storeC12];
                end
                    storeC12 = C12;
                end
                    storeC11=C11;
                end
                    storeC10=C10;
                end
                    storeC9=C9;
                end
                    storeC8=C8;
                end
                    storeC7=C7;
                end
                    storeC6=C6;
                end
                    storeC5=C5;
                end
                    storeC4=C4;
                end
                    storeC3=C3;
                end
                    storeC2=C2;
                end
                storeC1=C1;
                
                i=i+2;
            end
            if n == 1
                mrep=C1;
            end
            if n == 2
                mrep=C2;
            end
            if n == 3
                mrep=C3;
            end
            if n == 4
                mrep=C4;
            end
            if n == 5
                mrep=C5;
            end
            if n == 6
                mrep=C6;
            end
            if n == 7
                mrep=C7;
            end
            if n == 8
                mrep=C8;
            end
            if n == 9
                mrep=C9;
            end
            if n == 10
                mrep=C10;
            end
            if n == 11
                mrep=C11;
            end
            if n == 12
                mrep=C12;
            end
            if n == 13
                mrep=C13;
            end
            if n == 14
                mrep=C14;
            end
        end
        
        function out = plus(obj1,obj2)                           % Addition
            
            % MATLAB accepts putting multicomplex numbers inside arrays, so
            % the for loop inside the addition and substraction operator
            % overloading allows for the elementswise addition and
            % substraction of same sized arrays of multicomplex numbers.
            
            siz=size(obj1);
            %outfor = zeros(1,siz(2));
            for i=1:siz(2)
                [obj11,obj21] = consimulti(obj1(i),obj2(i));
                outfor(i)=multicomplex(obj11.zn+obj21.zn);
            end
            out=outfor;
        end
        
        function out = minus(obj1,obj2)                           % Substraction
            
            % Same algorithm structure as the addition operator.
            
            siz=size(obj1);
            %outfor = zeros(1,siz(2));
            for i=1:siz(2)
                [obj11,obj21] = consimulti(obj1(i),obj2(i));
                outfor(i) = multicomplex(obj11.zn-obj21.zn);
            end
            out=outfor;
        end
        
        function out = uminus(obj1)                               % Unary Minus
            out=multicomplex(-obj1.zn);
        end
        
        function out = uplus(obj1)                                % Unary Plus
            out=multicomplex(+obj1.zn);
        end
        
        function out = lt(self,other)              % Less than, self < other
            out = false;
            [self,other] = consimulti(self,other);
            if self.zn(1) < other.zn(1)
                out = true;
            end
        end
        
        function out = gt(self,other)           % Greater than, self > other
            out = false;
            [self,other] = consimulti(self,other);
            if self.zn(1) > other.zn(1)
                out = true;
            end
        end
        
        function out = le(self,other)    % Less than or equal, self <= other
            out = false;
            [self,other] = consimulti(self,other);
            if self.zn(1) <= other.zn(1)
                out = true;
            end
        end
        
        function out = ge(self,other) % Greater than or equal, self >= other
            out = false;
            [self,other] = consimulti(self,other);
            if self.zn(1) >= other.zn(1)
                out = true;
            end
        end
        
        function out = eq(self,other)              % Equality, self == other
            out = false;
            [self,other] = consimulti(self,other);
            if self.zn == other.zn
                out = true;
            end
        end
        
        function out = ne(self,other)             % Not equal, self ~= other
            out = true;
            [self,other] = consimulti(self,other);
            if self.zn == other.zn
                out = false;
            end
        end
        
        function out = mtimes(self,other)                         % Multiplication
            
            [self,other] = consimulti(self,other);
            %out=multicomplex(arrayM(matrep(self)*matrep(other))); % operation using built in matrices
            out = multicomplex(arr4matmulti(self,other)); % matrixless operation
        end
        
        function out = mrdivide(self,other)                       % Division
                      
            [self,other] = consimulti(self,other);
            
            % The following check for zero divisors is only turned on if
            % specified by the user
            if other.check_Z0 == 'T'
                if det(matrep(other)) == 0
                    error('Attempting to divide by a zero divizor multicomplex number')
                end
            end
            
            out=multicomplex(arrayM(matrep(self)/matrep(other)));
            
        end
        
        function out = mpower(input,p)                             % Power
            
            % First we check and see whether the power is itself a multicomplex
            % number, in which case the solution is of the form e^*log(input)
            
            l=length(input.zn);
            
            if isa(p,'multicomplex') == 1
    
                [input,p] = consimulti(input,p);
                
                out=exp(p*log(input));
                
            elseif imag(p) ~= 0
                
                input = log(input);
                
                input1 = cos(input * imag(p)) + sin(input * imag(p))*i;
                
                input2 = exp(input*real(p));
                
                out = input1*input2;
                
            elseif rem(p,1) == 0
                
                % if the power is an integer, the easiest method of multiplication is to
                % use the standard matrix multiplication p times.
                
                out=multicomplex(arrayM(matrep(input)^p)); %matrix
                %out=multicomplex(arr4mat(input,p)); %matrixless
            
            else
                
                % However, if it is a fractional power using an iterative de Moivres theorem is employed.
                % This algorithm breakes down the n order multicomplex number into the
                % two n-1 components and performs de Moivres. When you reach C_1, the
                % built in atan function is employed.
                if l == 1  % This just checks for a real number input
                    
                    out = multicomplex(input.zn^p);
                    
                else

                    out = multicomplex(arrayM(fractionalpow(matrep(input),p))); % for fractional power
                    %out = multicomplex(fracpow((input),p))
                    %out = multicomplex(cppfracpow_mex((input.zn),p));
                end
            end
        end
              
   %% Basic Function Overloading %%
        
        
        function out = atan2(input2,input1)                  % Arctan2
            
            % Works in a similar way as the other functions, breaking down
            % the C_n multicomplex number in to C_n-1 until you reach C_1.
            
            l=length(input1.zn);
            
            if l == 1   % Again checking for a real number input.
                
                out=multicomplex(atan2(input2.zn,input1.zn));
                
            else
                
                % The atan2 version of arctan is coded in here to take into
                % account the signs of the inputs.
                
                input=input2/input1;
                
                % Here we also have the approximation we derived for very
                % small imginary components.
                
                if  abs(input.zn((l/2)+1)) < 10^-4
                    
                    if input.check_convergence == 'T' && modcheck(multicomplex(input.zn(1:l/2)))
                        % Convergence criteria only turned on if
                        % specified by the user
                        
                        error('Convergence criteria for atan not met')
                    end
                    
                    half=l/2;
                    
                    outp1=atan2(multicomplex(input.zn(1:half)),multicomplex(1));
                    outp2=multicomplex(input.zn(half+1:end))/(1+(multicomplex(input.zn(1:half)))^2);
                    
                    out=multicomplex([outp1.zn(1:end),outp2.zn(1:end)]);
                    
                else
                    
                    input=((((input1^2)+(input2^2))^0.5)-input1)/input2;
                    
                    l=length(input.zn);
                    
                    ar=zeros(1,l);
                    ar((l/2)+1)=1;
                    ar=multicomplex(ar);
                    
                    outp=ar*0.5*log(multicomplex(arrayM((matrep(ar+input))/(matrep(ar-input)))));
                    out=2*outp;
                    
                end
                
            end
            
        end
        
        function out = atan(input)                           % Arctan
            
            l=length(input.zn);
            
            if l == 1
                
                out=multicomplex(atan(input.zn));
                
            elseif l==2
                
                
                outp=atan(input.zn(1)+(1i*input.zn(2)));
                out=multicomplex([real(outp),imag(outp)]);
                
            else
                
                
                if  abs(input.zn((l/2)+1)) < 10^-7
                    
                    half=l/2;
                    
                    outp1=atan2(multicomplex(input.zn(1:half)),multicomplex(1));
                    outp2=multicomplex(input.zn(half+1:end))/(1+(multicomplex(input.zn(1:half)))^2);
                    
                    out=multicomplex([outp1.zn(1:end),outp2.zn(1:end)]);
                    
                else
                    
                    ar=zeros(1,l);
                    ar((l/2)+1)=1;
                    ar=multicomplex(ar);
                    
                    out=ar*0.5*log(multicomplex(arrayM((matrep(ar+input))/(matrep(ar-input)))));
                    
                end
                
            end
        end
        
        function out = log(input)                               % Log
            
            l=length(input.zn);
            
            if l == 1
                
                out=multicomplex(log(input.zn));
                
            elseif l == 2
                
                % The correction below is for the cases where you need to calculate log(1+x), and x is very small
                % It uses the built in log1p(x) ? x which essentially eliminates
                % the roundoff error in the log(1+x) calculation.
                
                if abs(input.zn(1)) > 0.999999999 && abs(input.zn(1)) < 1.0000000001 && abs(input.zn(2)) < 10^-10
                    if input.zn(2) > 0
                        out=multicomplex([0.5*log1p(input.zn(2)^2),atan2(input.zn(2),input.zn(1))]);
                    else
                        out=multicomplex([-0.5*log1p(input.zn(2)^2),atan2(input.zn(2),input.zn(1))]);
                    end
                else
                    out=multicomplex([0.5*log(input.zn(1)^2+input.zn(2)^2),atan2(input.zn(2),input.zn(1))]);
                end
                
            else
                
                half=l/2;
                L_1 = multicomplex(input.zn(1:half));
                L_2 = multicomplex(input.zn(half+1:end));
                
                O_1=log(((L_1^2)+(L_2^2))^0.5);
                O_2=atan2(L_2,L_1);
                
                [O_1,O_2] = consimulti(O_1,O_2);
                out=multicomplex([O_1.zn,O_2.zn]);
                
            end
        end
        
        
        %       All of the trigonometric functions work similarly, separating the
        %       C_n multicomplex numbers into two C_n-1 coponents, extending the
        %       trigonometric definitions of C_1 complex numbers into the
        %       multicomplex domain.
        
        
        
        function out = sin(input)                              % Sin
            
            l=length(input.zn);
            
            if l == 1
                out=sin(input.zn);
                
                
            elseif l==2
                
                out=multicomplex([(sin(input.zn(1))*cosh(input.zn(2))),(cos(input.zn(1))*sinh(input.zn(2)))]);
                
            else
                
                half=l/2;
                L_1 = multicomplex(input.zn(1:half));
                L_2 = multicomplex(input.zn(half+1:end));
                
                O_1=sin(L_1)*cosh(L_2);
                O_2=cos(L_1)*sinh(L_2);
                
                out=multicomplex([O_1.zn,O_2.zn]);
                
            end
        end
        
        
        function out = cos(input)                              % Cos
            
            l=length(input.zn);
            
            if l == 1
                
                out=cos(input.zn);
                
            elseif l==2
                
                out=multicomplex([(cos(input.zn(1))*cosh(input.zn(2))),-(sin(input.zn(1))*sinh(input.zn(2)))]);
                
            else
                
                half=l/2;
                L_1 = multicomplex(input.zn(1:half));
                L_2 = multicomplex(input.zn(half+1:end));
                
                O_1=cos(L_1)*cosh(L_2);
                O_2=-sin(L_1)*sinh(L_2);
                
                out=multicomplex([O_1.zn,O_2.zn]);
                
            end
        end
        
        
        function out = sinh(input)                        % Sinh
            
            l=length(input.zn);
            
            if l == 1
                
                out=sinh(input.zn);
                
                
            elseif l==2
                
                out=multicomplex([(sinh(input.zn(1))*cos(input.zn(2))),(cosh(input.zn(1))*sin(input.zn(2)))]);
                
            else
                
                half=l/2;
                L_1 = multicomplex(input.zn(1:half));
                L_2 = multicomplex(input.zn(half+1:end));
                
                O_1=sinh(L_1)*cos(L_2);
                O_2=cosh(L_1)*sin(L_2);
                
                out=multicomplex([O_1.zn,O_2.zn]);
                
            end
        end
        
        
        function out = cosh(input)                             % Cosh
            
            l=length(input.zn);
            
            if l == 1
                
                out=cosh(input.zn);
                
                
            elseif l==2
                
                out=multicomplex([(cosh(input.zn(1))*cos(input.zn(2))),(sinh(input.zn(1))*sin(input.zn(2)))]);
                
            else
                
                half=l/2;
                L_1 = multicomplex(input.zn(1:half));
                L_2 = multicomplex(input.zn(half+1:end));
                
                O_1=cosh(L_1)*cos(L_2);
                O_2=sinh(L_1)*sin(L_2);
                
                out=multicomplex([O_1.zn,O_2.zn]);
            end
        end
        
        
        function out = tan(input)                             % Tan
            
            out=sin(input)/cos(input);
            
        end
        
        
        function out = asin(input)                            % Arcsin
            
            % Note that the value of the real part of the input is restricted
            % to |x|<1.
            
            l=length(input.zn);
            
            if l == 1
                
                out=asin(input.zn);
                
                
            elseif l==2
                
                outp=asin(input.zn(1)+(1i*input.zn(2)));
                out=multicomplex([real(outp),imag(outp)]);
                
            else
                
                
                if  abs(input.zn((l/2)+1)) < 10^-4
                    
                    half=l/2;
                    
                    outp1=asin(multicomplex(input.zn(1:half)));
                    outp2=multicomplex(input.zn(half+1:end))/((1-(multicomplex(input.zn(1:half)))^2)^0.5);
                    
                    out=multicomplex([outp1.zn(1:end),outp2.zn(1:end)]);
                    
                else
                    
                    half=l/2;
                    ar=zeros(1,l);
                    ar((l/2)+1)=1;
                    ar=multicomplex(ar);
                    d=1-(input^2);
                    
                    d1=multicomplex(d.zn(1:half));
                    d2=multicomplex(d.zn(half+1:end));
                    
                    dd1=(d1^2+d2^2)^0.5;
                    dd2=atan2(d2,d1);
                    
                    out=-ar*log((ar*input)+(((dd1)^0.5)*exp(0.5*ar*dd2)));
                    
                end
            end
        end
        
        
        function out = acos(input)                          % Arccos
            
            % Note that the value of the real part of the input is restricted
            % to |x|<1.
            
            l=length(input.zn);
            
            if l == 1
                
                out=acos(input.zn);
                
                
            elseif l==2
                
                outp=acos(input.zn(1)+(1i*input.zn(2)));
                out=multicomplex([real(outp),imag(outp)]);
                
            else
                
                if  abs(input.zn((l/2)+1)) < 10^-4
                    
                    out=(pi/2)-asin(input);
                    
                else
                    
                    half=l/2;
                    ar=zeros(1,l);
                    ar((l/2)+1)=1;
                    ar=multicomplex(ar);
                    d=1-(input^2);
                    
                    d1=multicomplex(d.zn(1:half));
                    d2=multicomplex(d.zn(half+1:end));
                    
                    dd1=(d1^2+d2^2)^0.5;
                    dd2=atan2(d2,d1);
                    
                    out=-ar*log(input+(ar*(((dd1)^0.5)*exp(0.5*ar*dd2))));
                    
                end
            end
        end
        
           
        function out = exp(input)                          % Exp
            
            l=length(input.zn);
            
            if l == 1
                
                out=exp(input.zn);
                
            elseif l==2
                
                outp=exp(input.zn(1)+(1i*input.zn(2)));
                out=multicomplex([real(outp),imag(outp)]);
                
            else
                
                half=l/2;
                L_1 = multicomplex(input.zn(1:half));
                L_2 = multicomplex(input.zn(half+1:end));
                
                O_1=exp(L_1)*cos(L_2);
                O_2=exp(L_1)*sin(L_2);
                
                out=multicomplex([O_1.zn,O_2.zn]);
                
                
            end
        end
        
        %% Utility functions %%
        function out = real(C)
            % This function extracts the real part of a multicomplex number or multicomplex matrix
            for i = 1:size(C,1)
                for j = 1:size(C,2)
                    out(i,j) = C(i,j).zn(1);
                end
            end                  
        end
        
        function out = imgn(C)
            
            % This function essentially extracts the 'last' imaginary term, or the
            % one that contains all of the imaginary terms.
            
            out=C.zn(end);
            
        end
        
        function out = CXN(C,imgn_parts)
            % This function extracts the coefficient of the user inputed imaginary parts contained within imgn_parts array.
            % Ex: CXn(C,[1 2 4]) will extract the i1i2i4 coefficient of C. As the imaginary parts are commutitative e.g. i1i2i4 =
            % i4i1i2, the ordering of the parts in imgn_parts array does not matter.

            if isnumeric(C) || max(imgn_parts) > log2(length(C.zn)) || any(imgn_parts < 1)
                error('input not in required form or out of bounds')
            else
                index = 1;
                for n = imgn_parts
                    index = index + 2^(n-1);
                end
                out = C.zn(index);
            end
        end
        
        function out = CX2(C,im,in) 
            % Function to enable legacy code using CX2 function work with new CXN function
            out = CXN(C, [im, in]);
        end
        
        function [self,other] = consimulti(self,other)
            
            % Verifies that self and other are multicomplex, converts different sized
            % multicomplex numbers into the same format. This funciton enables the
            % operations between differently sized multicomplex numbers (and with real
            % numbers).
            
            if length(self)~=1 || length(other)~=1
                error('inputs to function cannot be arrays or matrices')
            end
            
            if isa(self,'double')==1 && isa(other,'multicomplex')==1
                
                aa=size(other.zn)-1;
                self = multicomplex([self,zeros(1,aa(2))]);
                
            elseif isa(self,'multicomplex')==1 && isa(other,'double')==1
                
                bb=size(self.zn)-1;
                other = multicomplex([other,zeros(1,bb(2))]);
                
            elseif isa(self,'multicomplex')==1 && isa(other,'multicomplex')==1
                
                aa=size(other.zn);
                bb=size(self.zn);
                
                if aa(2) > bb(2)
                    
                    self = multicomplex([self.zn,zeros(1,(aa(2)-bb(2)))]);
                    
                elseif aa(2)<bb(2)
                    
                    other = multicomplex([other.zn,zeros(1,(bb(2)-aa(2)))]);
                    
                else
                    
                end
            else
                
                error('Inputs to function are not of double or multicomplex type objects')
                
            end
            
        end
        
        
        function out = modcheck(input)
            
            % Function to be used with atan2 function for the convergence test,
            % which is needed in order to use the small imaginary term
            % approximation. Due to the nature of the small imaginary component
            % multicomplex numbers, when employing atan2 the convergence criteria
            % is almost always satisfied, so this check should be turned off for
            % increase in computational efficiency.
            
            % modc(input) is our first type of modulus
            % modc2(input) is our second type of modulus
            if (modc(input) < 1) && (modc2(input^2)<modc2(input))
                % Only of both of these are true the multicomplex number is likely
                % to converge
                
                out=('converge');
                
            else
                
                out=('diverge');
                
            end
            
        end
        
        
        function out = modc(input)
            % First type of modulus: recursive
            l=length(input.zn);
            input1=input;
            
            % input1 will be our first type of modulus (recursive)
            while l > 2
                
                l=length(input1.zn);
                input1=(((multicomplex(input1.zn(1:l/2)))^(2))+(multicomplex(input1.zn(l/2+1:end)))^(2))^(1/2);
                
            end
            
            out=input1.zn;
            
        end
        
        
        function out = modc2(input)
            % Second type of modulus: square root of all the coefficients of the
            % multicomplex object squared.
            
            out=(sum((input.zn(1:end)).^(2))^(1/2));
            
            
        end
        
        
        function out = conj(input)
            % Complex conjugate
            
            half = length(input.zn)/2;
            
            out = multicomplex([input.zn(1:half),-input.zn(half+1:end)]);
            
        end
        
        function out = repr(self)
            % This function returns a string representation of the multicomplex number
                        
            imag_strings = strings(length(self.zn),1);
            complex_order = log(length(self.zn))/log(2);
            
            for n = 1:complex_order
                imag_strings(2^(n-1)+1) = "i" + n;
                for i = 2^(n-1)+2: 2^(n)
                    imag_strings(i) = string(imag_strings(i-2^(n-1))) + string(imag_strings(2^(n-1)+1));
                end
            end
            
            string_mat = strings(0);
            if self.zn(1) ~= 0
                string_mat(1) = string(self.zn(1));
            end
            for i = 2:length(self.zn)
                if self.zn(i) ~= 0
                    string_mat(end+1) = string(self.zn(i)) + " * " + string(imag_strings(i));                  
                end
            end
            
            if isempty(string_mat)
                out = "0";  
            else
                out = join(string_mat, "  +  ");
            end
        end
        
    end
    
end