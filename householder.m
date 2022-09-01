clear all;

n = 4;
A = [ 4,1,-2,2 ; 1,2,0,1 ; -2,0,3,-2 ; 2,1,-2,-1 ]; 

%n = 4;
%A = hilb(n);

% ======================================================================
% La première étape consiste à tridiagonaliser la matrice A à l'aide de la méthode de Householder. 
% Ensuite, on va itérer sur le processus qui permet d'obtenir une factorisation
% QR afin de déterminer les valeurs propres (en repartant de la matrice T tridiagonale, semblable à A).
% ======================================================================

disp("===== PHASE 1: Tridiagonalisation =====");
[Q1, T] = triangularization(A, n);
T = arrondir_a_zero(T);

% AFFICHER LES RESULTATS
disp(T);
disp("On peut verifier a l aide de la fonction eig (pour l'instant), que A et T sont semblables:");
disp("Valeurs propres de A:");
disp(transpose(eig(A)));
disp("Valeurs propres de T:");
disp(transpose(eig(T)));

% ======================================================================

% DECOMMENTER LA LIGNE SUIVANTE POUR CONSTATER QUE LA TRIDIAG EST CONTRE-PRODUCTIVE
% T = A;

disp("===== PHASE 2: Valeurs propres de A =====");
lambda = eigenvalues(T, n);
disp(lambda);

% ======================================================================

disp("===== PHASE 3: Vecteurs propres de A =====");
for i = 1:n
	V = eigenvectors(Q1, T, lambda(i), n);
	disp("Vecteur propre : ");
	disp(V);
end
   
% ======================================================================

function result = my_sign(x)
    %FONCTION SIGNE LEGEREMENT MODIFIEE AVEC : SIGN(0)=1
	result = sign(x);
	if (0 == result)
		result = 1;
    end
end

function [Q, T] = triangularization(A, n)
    % https://en.wikipedia.org/wiki/Householder_transformation#Tridiagonalization

	if (1 ~= issymmetric(A))
		error("A doit etre symétrique!\n");
    end

	I = eye(n); Q = eye(n); T = A;

	for c = 1:(n-1)
        % DECOMPOSITION DE T EN UN VECTEUR COLONNE U
		u = T(:,c);
		u(1:c) = 0; % A FAIRE EN PREMIER POUR NE PAS FAUSSER norm(u)
		alpha = - my_sign(u(c+1)) * norm(u);
		r = sqrt(0.5 * alpha * (alpha - u(c+1)));
		u(c+1) = (u(c+1) - alpha) / (2 * r);
		u(c+2:n) = u(c+2:n) / (2 * r);

		H = I - 2 * u * transpose(u) / norm(u)^2;

		Q = Q * H;
		T = transpose(H) * T * H;
	end
end

function [Q, R] = decompose_QR(A, n)
	I = eye(n); Q = eye(n); R = A;

	for c = 1:(n-1)
        % DECOMPOSITION DE R EN UN VECTEUR COLONNE U
		u = R(:, c);

		alpha = -sign(u(c));
		u( 1:(c-1) ) = 0;
		u(c) = u(c) - alpha * norm(u);

		H = I - 2 * u * transpose(u) / norm(u)^2;

		Q = Q * H;
		R = H * R;
    end

	R = transpose(Q) * A;
end

function [lambda] = eigenvalues(T, n)
    prev = zeros(n, 1);
	seuil = 1e-6;

	for i = 1:500
		[Q,R] = decompose_QR(T, n);
		T = R* Q;
        
        curr = diag(T);
		if (norm(prev-curr) < seuil)
			disp("Nombre d iterations necessaires");
            disp(i);
			break;
        end
		prev = curr;

	end
	lambda = transpose(diag(T));
end

function V = eigenvectors(Q, T, lambda, n)
    % https://core.ac.uk/download/pdf/82190245.pdf (7. Sturm sequences)

	V = ones(n, 1);
	prodB = 1; % PRODUIT DES COEFFICIENTS "b"
	P0 = 1;
	P1 = T(1,1) - lambda;

	for i = 2:n
		prodB = prodB * T(i,i-1);
		V(i) = power((-1),(i-1)) * P1 / prodB;
		P2 = (T(i,i) - lambda) * P1 - pow2(T(i,i-1)) * P0;
		P0 = P1;
		P1 = P2;
    end

    % PERMET D OBTENIR LES VECTEURS PROPRES DE A :
	V = Q * V;
end 

function result = arrondir_a_zero(A)
	% FORCER LES VALEURS PROCHES DE ZERO A ETRE NULLES
	seuil = 1e-12;
	result = A;
	result(abs(result) < seuil) = 0;
end
