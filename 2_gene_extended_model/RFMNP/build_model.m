function [Elong, Genes, Init, mRNA, errori] = build_model(codonPath, ...
    seqData, chunkSize, initRatePath, mRNAPath)
% [Elong, Genes, Init, mRNA, errori] = build_model(codonPath, ...
%    seqPath, chunkSize, initRatePath, mRNAPath)
%
%   Generating [Elong] (elongation rates per chunk) vector for genes based
%   on coding sequence, codon translation times & chunk size.
%
% Input:
%   [codonPath] file containing codon decoding times
%   [seqData] fasta file with gene sequences / table with vars (name, seq)
%   [chunkSize] size of RFM sites along the transcript
%   [initRatePath] (optional) file containing initiation rates
%   [mRNAPath] (optional) file containing mRNA levels
%
% Method:
%   Let N denote the size of a chunk in codons;
%   Go over the coding sequence and compute for each *non-overlapping*
%   chunk i of size N the total time Ti it takes to translate that chunk by
%   *summing* the times.
%   of all the codons in this chunk.  The rate (lambda) of the i-th chunk
%   is 1/Ti.
%   * Note that the last chunk may be shorter than N (if the length of the
%   coding sequence does not divide by N).
%   For example if the sequence is of length 12 codons and the cunk size is
%   5, there are 3 chunks: 1) codons 1-5; 2) codons 6-10; 3) codons 11-12.
%
% Hadas Zur / Tuller Lab
% (adapted from generateLamdot.m)
%
% Changelog:
%   21/02/16: parsing mRNA files. some additional input options,
%   e.g. providing sequences and codon times as variables. (Alon)
%
% *****************
%
% For example, if the input is the following file (note that a file can include many genes in fasta format)
% >zur1
% ATGAACAATATTTATATTATTATATATACTATCTCTATC
%
% With a time Tx for each coding codon x (need 61 times.. we ignore the STOP codons.. if there are stop codons in the sequence we should return 'error') 
% and the chunk size is 5, and the rate is l_init
% We perform the following steps:
% Compute the Lamdot:
% l_1 = 1/( T_ATG + T_AAC + T_AAT + T_ATT + T_TAT), l_2 = 1/(T_ATT + T_ATT + T_ATA + T_TAT + T_ACT), l_3 = 1/(T_ATC + T_TCT + T_ATC)
% Now run the RFM with lambda and lambdot to get the rate: R_z and the "densities" for each chunk (the vector of pi): p_z1, p_z2, p_z3

    global NC; NC = 61;
    global stopCodons; stopCodons = {'TAG'; 'TAA'; 'TGA'};

    Elong = [];
    Genes = [];
    Init = [];
    mRNA = [];
    errori = 0;

    if nargin < 5
        mRNAPath = false;
    end
    if nargin < 4
        initRatePath = false;
    end
    
%     fid = fopen('RFMLog.txt', 'a');

    % parse codon translation times
    [codons, ctimes, success] = parseCodons(codonPath);
    if (~success)
        errori = 1;
        return;
    end

    if ischar(seqData)
        % parse coding sequences file
        [success, genenames, seqs] = parseFasta(seqData);
        if (~success)
		    errori = 2;
            return;
        end
        % (keeping legacy variables and generating a new table)
        Genes = RNASeq(seqs, genenames);
    else
        % coding sequences are already in RNASeq array [seqData]
        genenames = {seqData.Name};
		for i=1:length(seqData)
			seqData(i).Bases = strcat(seqData(i).Codons{:});
		end
        seqs = {seqData.Bases}; %Whenever I save the RNASeq I remove the bases to save space so this no longer works on save data without the above for loop
        Genes = copy(seqData);

    end
    
    % parse initiation rate file
    if initRatePath
        [Init, success] = parseInitRate(initRatePath, genenames);
        if (~success)
            errori = 3;
            return;
        end
    end
    
    % parse mRNA file
    if mRNAPath
        [mRNA, success] = parseInitRate(mRNAPath, genenames);
        if (~success)
            errori = 4;
            return;
        end
    end
    
    Elong = cell(length(genenames), 1);
    for i = 1:length(genenames)
        [Elong{i}] = segmentChunks(seqs{i}, codons, ctimes, chunkSize);
    end

end

%%

function [lamdot] = segmentChunks(seq, codons, ctimes, chunkSize)

    % lamdot = -1: the last codon is not a Stop codon
    % lamdot = -2: the coding sequence contains a stop codon
    % lamdot = -3: the coding sequence contains an illegal codon
    % lamdot = -4: the coding sequence is shorter than the chunk size
    % lamdot = -5: the coding sequence is not divisible by 3 and hence not composed of legal codons
    % lamdot = -6: unkown error occured

    global stopCodons;
    
    try
    
        if (isempty(find(strcmp(seq(end-2:end), stopCodons), 1)))
            lamdot = -1;
            return;
        end
        seq = seq(1:end-3);

        if ((length(seq)/3)/chunkSize < 1)
            lamdot = -4;
            return;
        end
        if (rem(length(seq), 3))
            lamdot = -5;
            return;
        end

        R = rem(length(seq)/3, chunkSize)/chunkSize;
        numChunks = floor((length(seq)/3)/chunkSize);
        if (numChunks == 0)
            numChunks = 1;
        end
        chunkSize = chunkSize*3; % convert to nts
        M = 0;
        if (R ~= 0)
            M = 1; % an extra, partial chunk
        end
        lamdot = nan(numChunks+M,1);
        count = 1;
        for (i=1:chunkSize:((numChunks*chunkSize)+M))
            lamda = 0;
            % get chunk sequence
            if (i == ((numChunks*chunkSize)+M))
                % last chunk
                if (R ~= 0)
                    % partial chunk
                    cseq = seq(i:i+R*chunkSize-1);
                else
                    fprintf('test');
                end
            else
                cseq = seq(i:i+chunkSize-1);
            end
            for (j=1:3:length(cseq)) % for each codon
                ind = find(strcmp(cseq(j:j+2), codons));
                if (~isempty(ind))
                    lamda = lamda + ctimes(ind);
                elseif (~isempty(find(strcmp(cseq(j:j+2), stopCodons), 1)))
                    lamdot = -2;
                    return;
                elseif (isempty(find(strcmp(cseq(j:j+2), stopCodons), 1)))
                    lamdot = -3;
                    return;
                end
            end
            lamdot(count) = 1/lamda;
            count = count+1;
        end
    catch err
        lamdot = -6;
        return;
    end
end

%%

function [codons, ctimes, success] = parseCodons(file)

    global NC;
    
    success = 1;
    codons = [];
    ctimes = [];
    
%     try
        fid = fopen(file);    
        C = textscan(fid, '%s\t%s');
        fclose(fid);

        codons = C{1,1};
        ctimes = C{1,2};

        if ((length(codons) == NC) && (length(ctimes) == NC))
            [success, ctimes] = validateCodons(codons, ctimes);
            ctimes = cellfun(@str2num, ctimes);
        end
%     catch err
%         success = 0;
%     end
    
end

%%

function [success, ctimes] = validateCodons(codons, ctimes)

    global NC;
    global stopCodons
    
    vcodons = setxor(stopCodons, fieldnames(codoncount('')));
    
    success = 0;

    if length(codons) ~= length(unique(codons))
        return;
    elseif nnz(ismember(codons, vcodons)) ~= NC
        return;
    else
        success = 1;
    end
    
end

%%

function [success, genes, seqs] = parseFasta(file)

    genes = {};
    seqs = {};
    
    try
        fid = fopen(file);
        l = textscan(fid, '%s', 1, 'Delimiter', '\n');
        seq = [];
        while (~cellfun('isempty',l))

            l = l{1}{1};

            if (strncmp('>', l, 1))
                if (~isempty(seq))
                   if (~strcmp('Sequence unavailable', seq))
                        genes{end+1,1} = gene;
                        seqs{end+1,1} = seq;               
                   end
                   seq = [];
                end
                gene = strtrim(l(2:end));
                if (isempty(gene))
                    success = 0;
                    return;
                end
                
                l = textscan(fid, '%s', 1, 'Delimiter', '\n');
            else
                    seq = [seq, l];
                    l = textscan(fid, '%s', 1, 'Delimiter', '\n');
            end
        end
        if (~isempty(seq)) % check if we ignored anything (check: does l contain anything?)
            genes{end+1,1} = gene;
            seqs{end+1,1} = seq;               
            seq = [];
        end

        fclose(fid);
        success = 1;
    catch err
		disp(err)
        success = 0;
    end

end

%%

function [initRates, success] = parseInitRate(file, genes)

    success = 0;
    
    initRates = nan(length(genes), 1);
    try
        fid = fopen(file);    
        C = textscan(fid, '%s\t%s');
        fclose(fid);

        rgenes = C{1,1};
        rates = C{1,2};
        
        if (isempty(rates))
            return;
        end
        
        if (length(rates) == 1)
            initRates = str2num(rates{1});
            success = 1;
            return;
        end
%         if ((length(rgenes) == length(genes)) && (length(rates) == length(genes)))
            for (i=1:length(rgenes))
                ind = find(strcmp(rgenes{i}, genes));
                if (~isempty(ind))
                    initRates(ind) = str2num(rates{i}); 
%                 else
%                     success = 0;
%                     return;
                end
            end
%             initRates = cellfun(@str2num, initRates);
            success = 1;
%         end
    catch err
        success = 0;
    end

end

%%

function [times] = fixStopTimes(codons, times)

    global stopCodons;
    
    for (i=1:length(stopCodons))
        idx = find(strcmp(stopCodons{i}, codons));
        times(idx) = 0;
    end
end
