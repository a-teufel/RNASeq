classdef RNASeq < dynamicprops & matlab.mixin.Copyable
%Basically this is so we can just instantiate an RNASeq object in the
%script and keep things neat

	properties (Access = public)
		Name      = '';
		Bases     = '';
		Density   = [];
		Elong     = [];
		Rate      = [];
		Pool      = [];
		Codons    = {''};
		AA        = '';
		score     = 1;
	end

	methods

		function obj = RNASeq(bases, name) %constructor. Adds the bases and cuts them into codons

			if iscell(bases)
				if nargin < 2
					name = cell(size(bases));
				end
				% create an array of RNASeq for each given sequence
				obj = cellfun(@RNASeq, bases, name, 'UniformOutput', false);
				obj = cat(1, obj{:});
				return
			end
			if nargin > 1
				obj.Name = name;
			end

			if ischar(bases) == 1
				obj.Bases = bases;
			else
				error('Value must be string')
			end
			CutSeq;
			obj.AA = nt2aa(bases);

		function CutSeq %cuts bases into codons
			i = 1;
			n = 1;
			while i < (length(obj.Bases))
				obj.Codons{1,n} = obj.Bases(i:(i+2));
				i = i+3;
				n = n+1;
			end
		end

		end

	function MutSeq = mutate(RNA1)
	%takes an object of the RNASeq class and alters its 'Codons' by changing the first indexed codon to a random location
		MutSeq = copy(RNA1);
		cod = randi(ceil(length(MutSeq.Codons)));  % select a codon uniformly
		while any(strcmpi(MutSeq.Codons{cod}, {'ATG', 'TGG', 'TGA'})) %Dont select the codons that don't have synonymous sibling codons
			cod = randi(length(MutSeq.Codons));
		end
		temp = ReplaceCodon(MutSeq, cod);
		MutSeq.Codons{cod} = temp;
		MutSeq.Bases((cod-1)*3+1 : cod*3) = temp;  % we keep bases updated as well
	end

	function NewCodon = ReplaceCodon(MutSeq, cod)
	% new function may replace not only the 3rd position of the codon,
	% so that AA with 6 degrees of freedom are covered
		OldCodon = MutSeq.Codons{cod};
		codonmap = codonbias('');
		synonyms = codonmap.(aminolookup(nt2aa(OldCodon))).Codon;
		synonyms(strcmp(OldCodon, synonyms)) = [];  % remove self
		if isempty(synonyms)
			error('no synonyms');
		end
		NewCodon = synonyms{randi(length(synonyms))};
	end

	function new_obj = new_prop(new_obj,name,value)
		if ischar(name)
			new_obj.addprop(name)
			new_obj.(name) = value;
		else
			disp('not a char')
		end
	end

	function MaxRate = CalcMaxRate(codonfile, Genes, ChunkSize, Init, mRNA, PoolSize, InitFactor)
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%Calculate our theoretical maximum protein production
	%We assume that the max pp rate is just all fast codons.
		%Replace with all fast codons
		for g =1:length(Genes)
			Genes(g).Codons(strcmp(Genes(g).Codons, 'TAC'))={'TAT'};
		end
		%Calculate PP Rate Using Alons Code
		ElongM = build_model(codonfile, Genes, ChunkSize);
		[~, Rate, ~, ~, ~] = solve_RFMNP(ElongM, Init, mRNA, PoolSize, InitFactor);
		MaxRate = sum(cellfun(@(x) x(end), Rate));
	end

	function [Genes, OrgScore] = ScoreFn(Genes, Rate, MaxRate)
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%Maximize Protein Production Rate while ensuring each gene
	%has the same portion of the protein production rate.
	%Effectively, PP_Rate(1) = PP_Rate(2) if we crank the W
		W = 0;  % weight of score A vs B
		weight = W/(1+W);
		TotalRate = sum(cellfun(@(x) x(end), Rate));  % current
		F1 = (TotalRate)/MaxRate; % Percentage of maximum score
		F2 = ((Rate{1}(end) - Rate{2}(end)).^2) / TotalRate.^2; %Expression Error (vs 0.5)
		OrgScore = F1*(1-weight) + weight*(1-F2); %Should be bounded between 0,1
	end

	function write_csv(obj,csv_file,varargin)
	%Transpose every property for csv writing onto a new struct
		temp_obj = copy(obj);
		temp_obj.Density  = obj.Density.';
		temp_obj.Elong    = obj.Elong.';
		temp_obj.Pool     = obj.Pool.';
		temp_obj.Rate     = obj.Rate.';
		temp_obj.AA       = '';
		temp_obj.Bases    = ''; %This should be empty to save space
		temp_obj.Codons   = {''}; %This should be empty to save space
		temp_obj.AA       = ''; %This should be empty to save space
		%Transpose any other property names that the user must pass in
		for arg = 1:length(varargin)
			temp_obj = temp_obj.new_prop(inputname(arg+2),varargin{arg}.');
		end
		struct2csv(temp_obj,csv_file)
		clear temp_obj
	end
end

	methods (Static)
	function Fit = FitFn(cur, mut,mRNA) %Probability that this mutation will be accepted. Dependant on protein output and expression error.
		if mut > cur
			Fit = 1;
		else
			Fit = exp((-2*500)*(cur-mut));
		end
	end


	function save_fig(image_file,y_data,x_label,y_label,title_text)
	% Saves a figure with x, y axis labels and a title
		fig = figure;
		plot(y_data)
		xlabel(x_label)
		ylabel(y_label)
		title(title_text)
		savefig(fig,image_file)
	end


	end
end
