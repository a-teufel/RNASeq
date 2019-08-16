function [ ramp_ratio, body_ratio ] = slow_codon_fraction( gene,codons )
    %Computes fraction of slow codons for the ramp and rest of
    %sequence
    ramp_num_codons = zeros(1,length(codons));
    body_num_codons = zeros(1,length(codons));
    for i=1:length(codons)
        for p=1:length(gene.Codons)
            if strcmp(gene.Codons{p},codons{i}) && p < length(gene.Codons)/10
                ramp_num_codons(i) = ramp_num_codons(i)+1;
		    end
            if strcmp(gene.Codons{p},codons{i})
				sprintf('Codon number: %i, Codon: %s',i,codons{i})
                body_num_codons(i) = body_num_codons(i)+1;
            end
        end
    end

    ramp_ratio = ramp_num_codons(1)/(ramp_num_codons(1)+ramp_num_codons(2));
    body_ratio = body_num_codons(1)/(body_num_codons(1)+body_num_codons(2));
end

