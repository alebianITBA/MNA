function melcoefficient = mel_mapping(value)
	melcoefficient = 1125 * log(1+(value/700));
end