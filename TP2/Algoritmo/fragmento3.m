function answer = fragmento3()

	
	input_user = parser();

	if (length(input_user)==0)
		printf('Error en los parametros. Revisar input.txt\n');
		return;
	endif

	
	

	printf('*************************************************************************************\n');
	printf('    Bienvenido, los Parametros con los que correra el reconocedor de voz \n');
	printf('    se pueden cambiar en el archivo input.txt, todas las pruebas y entrenamientos \n');
	printf('    correran en principio con los mismos parametros.  \n');
	printf('*************************************************************************************\n');

	%VARIABLES SECTION
		%pre-emphasize coefficient
		a = input_user{1} 
		%sample rate in hz
		fs = input_user{2} 
		%number of bits per sample
		bps = input_user{3} 
		%frame length in ms
		frame_length_ms = input_user{4} 
		%coefficient of overlap between 2 consecutive frames
		overlapping_coef = input_user{5}
		
		%Lower frecuency used in Mel filterBank
		lower_bank_filter_frec = input_user{6}
		%Higher frecuency used in Mel filterBank
		higher_bank_filter_frec = input_user{7}
		%160 samples as it was recorded at 8000hz with a 
		%frame_length_ms = 0,02s = 20ms to guarantee stationary
		frame_length = fs * (frame_length_ms/1000)
		%80 samples as it was recorded at 8000hz with a 
		%frame_length_ms = 0,01s = 10ms of overlap
		frame_length_overlaped = fs * overlapping_coef * (frame_length_ms/1000)
		%quantity of filter banks
		quantity_of_banks = input_user{8}

		nfft = input_user{9}
warning ('off', 'Octave:broadcast');

%TRAINING - FRAGMENT 1 - MEL COEFFICIENTS + VECTOR QUANTIZATION
[y,fs,bps] = wavread('../Entrenamiento/daniel_fragmento_1.wav');
mel_coef_entrenamiento_1 = mfcc(y,bps);
vector_quantization_1 = vq(mel_coef_entrenamiento_1,16);

[y,fs,bps] = wavread('../Entrenamiento/erika_fragmento_1.wav');
mel_coef_entrenamiento_2 = mfcc(y,bps);
vector_quantization_2 = vq(mel_coef_entrenamiento_2,16);

[y,fs,bps] = wavread('../Entrenamiento/eduardo_fragmento_1.wav');
mel_coef_entrenamiento_3 = mfcc(y,bps);
vector_quantization_3 = vq(mel_coef_entrenamiento_3,16);

[y,fs,bps] = wavread('../Entrenamiento/KevinK_fragmento_1.wav');
mel_coef_entrenamiento_4 = mfcc(y,bps);
vector_quantization_4 = vq(mel_coef_entrenamiento_4,16);

[y,fs,bps] = wavread('../Entrenamiento/helga_fragmento_1.wav');
mel_coef_entrenamiento_5 = mfcc(y,bps);
vector_quantization_5 = vq(mel_coef_entrenamiento_5,16);

[y,fs,bps] = wavread('../Entrenamiento/kevinj_fragmento_1.wav');
mel_coef_entrenamiento_6 = mfcc(y,bps);
vector_quantization_6 = vq(mel_coef_entrenamiento_6,16);

[y,fs,bps] = wavread('../Entrenamiento/nicolas_fragmento_1.wav');
mel_coef_entrenamiento_7 = mfcc(y,bps);
vector_quantization_7 = vq(mel_coef_entrenamiento_7,16);

[y,fs,bps] = wavread('../Entrenamiento/paul_fragmento_1.wav');
mel_coef_entrenamiento_8 = mfcc(y,bps);
vector_quantization_8 = vq(mel_coef_entrenamiento_8,16);

[y,fs,bps] = wavread('../Entrenamiento/jorge_fragmento_1.wav');
mel_coef_entrenamiento_9 = mfcc(y,bps);
vector_quantization_9 = vq(mel_coef_entrenamiento_9,16);

[y,fs,bps] = wavread('../Entrenamiento/alex_fragmento_1.wav');
mel_coef_entrenamiento_10 = mfcc(y,bps);
vector_quantization_10 = vq(mel_coef_entrenamiento_10,16);

name_vector = ['Daniel';'Erika';'Eduardo';'KevinK';'Helga';'KevinJ';'Nicolas';'Paul';'Jorge';'Alex'];
correct = 0;
quantity_of_tests=0;

printf('Iniciamos las pruebas...\n\n');
%%**************************************************************
%%**************************************************************

printf('Utilizando el Fragmento 3 para la prueba...\n');
% %FRAGMENT 3
[y,fs,bps] = wavread('../Pruebas/daniel_fragmento_3.wav');
mel_coef_prueba_daniel_3 = mfcc(y,bps);

[y,fs,bps] = wavread('../Pruebas/erika_fragmento_3.wav');
mel_coef_prueba_erika_3 = mfcc(y,bps);

[y,fs,bps] = wavread('../Pruebas/eduardo_fragmento_3.wav');
mel_coef_prueba_eduardo_3 = mfcc(y,bps);

[y,fs,bps] = wavread('../Pruebas/KevinK_fragmento_3.wav');
mel_coef_prueba_KevinK_3 = mfcc(y,bps);

[y,fs,bps] = wavread('../Pruebas/helga_fragmento_3.wav');
mel_coef_prueba_helga_3 = mfcc(y,bps);

[y,fs,bps] = wavread('../Pruebas/kevinj_fragmento_3.wav');
mel_coef_prueba_kevinj_3 = mfcc(y,bps);

[y,fs,bps] = wavread('../Pruebas/nicolas_fragmento_3.wav');
mel_coef_prueba_nicolas_3 = mfcc(y,bps);

[y,fs,bps] = wavread('../Pruebas/paul_fragmento_3.wav');
mel_coef_prueba_paul_3 = mfcc(y,bps);

[y,fs,bps] = wavread('../Pruebas/jorge_fragmento_3.wav');
mel_coef_prueba_jorge_3 = mfcc(y,bps);

[y,fs,bps] = wavread('../Pruebas/alex_fragmento_3.wav');
mel_coef_prueba_alex_3 = mfcc(y,bps);


%Prueba daniel
mean_distance(1) = meandist(mel_coef_prueba_daniel_3, vector_quantization_1);
mean_distance(2) = meandist(mel_coef_prueba_erika_3, vector_quantization_1);
mean_distance(3) = meandist(mel_coef_prueba_eduardo_3, vector_quantization_1);
mean_distance(4) = meandist(mel_coef_prueba_KevinK_3, vector_quantization_1);
mean_distance(5) = meandist(mel_coef_prueba_helga_3, vector_quantization_1);
mean_distance(6) = meandist(mel_coef_prueba_kevinj_3, vector_quantization_1);
mean_distance(7) = meandist(mel_coef_prueba_nicolas_3, vector_quantization_1);
mean_distance(8) = meandist(mel_coef_prueba_paul_3, vector_quantization_1);
mean_distance(9) = meandist(mel_coef_prueba_jorge_3, vector_quantization_1);
mean_distance(10) = meandist(mel_coef_prueba_alex_3, vector_quantization_1);

[min_value, index_of_min_value] = min(mean_distance);
if(strcmp('Daniel ', name_vector(index_of_min_value,:)))
	correct++;
endif
quantity_of_tests++;
printf('El individuo %s segun el reconocedor de voces es %s \n',' Daniel ', name_vector(index_of_min_value,:));

%Prueba erika
mean_distance(1) = meandist(mel_coef_prueba_daniel_3,vector_quantization_2);
mean_distance(2) = meandist(mel_coef_prueba_erika_3,vector_quantization_2);
mean_distance(3) = meandist(mel_coef_prueba_eduardo_3,vector_quantization_2);
mean_distance(4) = meandist(mel_coef_prueba_KevinK_3,vector_quantization_2);
mean_distance(5) = meandist(mel_coef_prueba_helga_3,vector_quantization_2);
mean_distance(6) = meandist(mel_coef_prueba_kevinj_3,vector_quantization_2);
mean_distance(7) = meandist(mel_coef_prueba_nicolas_3, vector_quantization_2);
mean_distance(8) = meandist(mel_coef_prueba_paul_3, vector_quantization_2);
mean_distance(9) = meandist(mel_coef_prueba_jorge_3, vector_quantization_2);
mean_distance(10) = meandist(mel_coef_prueba_alex_3, vector_quantization_2);

[min_value, index_of_min_value] = min(mean_distance);
if(strcmp('Erika  ', name_vector(index_of_min_value,:)))
	correct++;
endif
quantity_of_tests++;
printf('El individuo %s segun el reconocedor de voces es %s \n',' Erika ', name_vector(index_of_min_value,:));

%Prueba eduardo
mean_distance(1) = meandist(mel_coef_prueba_daniel_3,vector_quantization_3);
mean_distance(2) = meandist(mel_coef_prueba_erika_3,vector_quantization_3);
mean_distance(3) = meandist(mel_coef_prueba_eduardo_3,vector_quantization_3);
mean_distance(4) = meandist(mel_coef_prueba_KevinK_3,vector_quantization_3);
mean_distance(5) = meandist(mel_coef_prueba_helga_3,vector_quantization_3);
mean_distance(6) = meandist(mel_coef_prueba_kevinj_3,vector_quantization_3);
mean_distance(7) = meandist(mel_coef_prueba_nicolas_3, vector_quantization_3);
mean_distance(8) = meandist(mel_coef_prueba_paul_3, vector_quantization_3);
mean_distance(9) = meandist(mel_coef_prueba_jorge_3, vector_quantization_3);
mean_distance(10) = meandist(mel_coef_prueba_alex_3, vector_quantization_3);
[min_value, index_of_min_value] = min(mean_distance);
if(strcmp('Eduardo', name_vector(index_of_min_value,:)))
	correct++;
endif
quantity_of_tests++;
printf('El individuo %s segun el reconocedor de voces es %s \n',' Eduardo ', name_vector(index_of_min_value,:));

%Prueba KevinK
mean_distance(1) = meandist(mel_coef_prueba_daniel_3,vector_quantization_4);
mean_distance(2) = meandist(mel_coef_prueba_erika_3,vector_quantization_4);
mean_distance(3) = meandist(mel_coef_prueba_eduardo_3,vector_quantization_4);
mean_distance(4) = meandist(mel_coef_prueba_KevinK_3,vector_quantization_4);
mean_distance(5) = meandist(mel_coef_prueba_helga_3,vector_quantization_4);
mean_distance(6) = meandist(mel_coef_prueba_kevinj_3,vector_quantization_4);
mean_distance(7) = meandist(mel_coef_prueba_nicolas_3, vector_quantization_4);
mean_distance(8) = meandist(mel_coef_prueba_paul_3, vector_quantization_4);
mean_distance(9) = meandist(mel_coef_prueba_jorge_3, vector_quantization_4);
mean_distance(10) = meandist(mel_coef_prueba_alex_3, vector_quantization_4);

[min_value, index_of_min_value] = min(mean_distance);
if(strcmp('KevinK ', name_vector(index_of_min_value,:)))
	correct++;
endif
quantity_of_tests++;
printf('El individuo %s segun el reconocedor de voces es %s \n',' KevinK ', name_vector(index_of_min_value,:));

%Prueba Helga
mean_distance(1) = meandist(mel_coef_prueba_daniel_3,vector_quantization_5);
mean_distance(2) = meandist(mel_coef_prueba_erika_3,vector_quantization_5);
mean_distance(3) = meandist(mel_coef_prueba_eduardo_3,vector_quantization_5);
mean_distance(4) = meandist(mel_coef_prueba_KevinK_3,vector_quantization_5);
mean_distance(5) = meandist(mel_coef_prueba_helga_3,vector_quantization_5);
mean_distance(6) = meandist(mel_coef_prueba_kevinj_3,vector_quantization_5);
mean_distance(7) = meandist(mel_coef_prueba_nicolas_3, vector_quantization_5);
mean_distance(8) = meandist(mel_coef_prueba_paul_3, vector_quantization_5);
mean_distance(9) = meandist(mel_coef_prueba_jorge_3, vector_quantization_5);
mean_distance(10) = meandist(mel_coef_prueba_alex_3, vector_quantization_5);

[min_value, index_of_min_value] = min(mean_distance);
if(strcmp('Helga  ', name_vector(index_of_min_value,:)))
	correct++;
endif
quantity_of_tests++;
printf('El individuo %s segun el reconocedor de voces es %s \n',' Helga ', name_vector(index_of_min_value,:));

%Prueba KevinJHanna
mean_distance(1) = meandist(mel_coef_prueba_daniel_3,vector_quantization_6);
mean_distance(2) = meandist(mel_coef_prueba_erika_3,vector_quantization_6);
mean_distance(3) = meandist(mel_coef_prueba_eduardo_3,vector_quantization_6);
mean_distance(4) = meandist(mel_coef_prueba_KevinK_3,vector_quantization_6);
mean_distance(5) = meandist(mel_coef_prueba_helga_3,vector_quantization_6);
mean_distance(6) = meandist(mel_coef_prueba_kevinj_3,vector_quantization_6);
mean_distance(7) = meandist(mel_coef_prueba_nicolas_3, vector_quantization_6);
mean_distance(8) = meandist(mel_coef_prueba_paul_3, vector_quantization_6);
mean_distance(9) = meandist(mel_coef_prueba_jorge_3, vector_quantization_6);
mean_distance(10) = meandist(mel_coef_prueba_alex_3, vector_quantization_6);

[min_value, index_of_min_value] = min(mean_distance);
if(strcmp('KevinJ ', name_vector(index_of_min_value,:)))
	correct++;
endif
% subplot(3,3,1);
% plot(mean_distance,'*');
% title('Dist. minima-Entrenamiento: kevinj');
% xlabel('name_vector = ["Daniel"=1;"Erika"=2;"Eduardo"=3;"KevinK"=4;"Helga"=5];');
% ylabel('Mean Distance');
quantity_of_tests++;
printf('El individuo %s segun el reconocedor de voces es %s \n',' KevinJ ', name_vector(index_of_min_value,:));

%Prueba nicolas
mean_distance(1) = meandist(mel_coef_prueba_daniel_3,vector_quantization_7);
mean_distance(2) = meandist(mel_coef_prueba_erika_3,vector_quantization_7);
mean_distance(3) = meandist(mel_coef_prueba_eduardo_3,vector_quantization_7);
mean_distance(4) = meandist(mel_coef_prueba_KevinK_3,vector_quantization_7);
mean_distance(5) = meandist(mel_coef_prueba_helga_3,vector_quantization_7);
mean_distance(6) = meandist(mel_coef_prueba_kevinj_3,vector_quantization_7);
mean_distance(7) = meandist(mel_coef_prueba_nicolas_3, vector_quantization_7);
mean_distance(8) = meandist(mel_coef_prueba_paul_3, vector_quantization_7);
mean_distance(9) = meandist(mel_coef_prueba_jorge_3, vector_quantization_7);
mean_distance(10) = meandist(mel_coef_prueba_alex_3, vector_quantization_7);

[min_value, index_of_min_value] = min(mean_distance);
if(strcmp('Nicolas', name_vector(index_of_min_value,:)))
	correct++;
endif
% subplot(3,3,3);
% plot(mean_distance,'*');
% title('Dist. minima-Entrenamiento: nicolas');
% xlabel('name_vector = ["Daniel"=1;"Erika"=2;"Eduardo"=3;"KevinK"=4;"Helga"=5];');
% ylabel('Mean Distance');
quantity_of_tests++;
printf('El individuo %s segun el reconocedor de voces es %s \n',' Nicolas ', name_vector(index_of_min_value,:));

% %Prueba paul
mean_distance(1) = meandist(mel_coef_prueba_daniel_3,vector_quantization_8);
mean_distance(2) = meandist(mel_coef_prueba_erika_3,vector_quantization_8);
mean_distance(3) = meandist(mel_coef_prueba_eduardo_3,vector_quantization_8);
mean_distance(4) = meandist(mel_coef_prueba_KevinK_3,vector_quantization_8);
mean_distance(5) = meandist(mel_coef_prueba_helga_3,vector_quantization_8);
mean_distance(6) = meandist(mel_coef_prueba_kevinj_3,vector_quantization_8);
mean_distance(7) = meandist(mel_coef_prueba_nicolas_3, vector_quantization_8);
mean_distance(8) = meandist(mel_coef_prueba_paul_3, vector_quantization_8);
mean_distance(9) = meandist(mel_coef_prueba_jorge_3, vector_quantization_8);
mean_distance(10) = meandist(mel_coef_prueba_alex_3, vector_quantization_8);

[min_value, index_of_min_value] = min(mean_distance);
if(strcmp('Paul   ', name_vector(index_of_min_value,:)))
	correct++;
endif
% subplot(3,3,5);
% plot(mean_distance,'*');
% title('Dist. minima-Entrenamiento: paul');
% xlabel('name_vector = ["Daniel"=1;"Erika"=2;"Eduardo"=3;"KevinK"=4;"Helga"=5];');
% ylabel('Mean Distance');
quantity_of_tests++;
printf('El individuo %s segun el reconocedor de voces es %s \n',' Paul ', name_vector(index_of_min_value,:));


% %Prueba jorge
mean_distance(1) = meandist(mel_coef_prueba_daniel_3,vector_quantization_9);
mean_distance(2) = meandist(mel_coef_prueba_erika_3,vector_quantization_9);
mean_distance(3) = meandist(mel_coef_prueba_eduardo_3,vector_quantization_9);
mean_distance(4) = meandist(mel_coef_prueba_KevinK_3,vector_quantization_9);
mean_distance(5) = meandist(mel_coef_prueba_helga_3,vector_quantization_9);
mean_distance(6) = meandist(mel_coef_prueba_kevinj_3,vector_quantization_9);
mean_distance(7) = meandist(mel_coef_prueba_nicolas_3, vector_quantization_9);
mean_distance(8) = meandist(mel_coef_prueba_paul_3, vector_quantization_9);
mean_distance(9) = meandist(mel_coef_prueba_jorge_3, vector_quantization_9);
mean_distance(10) = meandist(mel_coef_prueba_alex_3, vector_quantization_9);

[min_value, index_of_min_value] = min(mean_distance);
if(strcmp('Jorge  ', name_vector(index_of_min_value,:)))
	correct++;
endif
quantity_of_tests++;
printf('El individuo %s segun el reconocedor de voces es %s \n',' Jorge ', name_vector(index_of_min_value,:));


% %Prueba Alex
mean_distance(1) = meandist(mel_coef_prueba_daniel_3,vector_quantization_10);
mean_distance(2) = meandist(mel_coef_prueba_erika_3,vector_quantization_10);
mean_distance(3) = meandist(mel_coef_prueba_eduardo_3,vector_quantization_10);
mean_distance(4) = meandist(mel_coef_prueba_KevinK_3,vector_quantization_10);
mean_distance(5) = meandist(mel_coef_prueba_helga_3,vector_quantization_10);
mean_distance(6) = meandist(mel_coef_prueba_kevinj_3,vector_quantization_10);
mean_distance(7) = meandist(mel_coef_prueba_nicolas_3, vector_quantization_10);
mean_distance(8) = meandist(mel_coef_prueba_paul_3, vector_quantization_10);
mean_distance(9) = meandist(mel_coef_prueba_jorge_3, vector_quantization_10);
mean_distance(10) = meandist(mel_coef_prueba_alex_3, vector_quantization_10);

[min_value, index_of_min_value] = min(mean_distance);
if(strcmp('Alex   ', name_vector(index_of_min_value,:)))
	correct++;
endif
quantity_of_tests++;
printf('El individuo %s segun el reconocedor de voces es %s \n',' Alex ', name_vector(index_of_min_value,:));

printf('La cantidad de matches son %d de %d\n\n',correct,quantity_of_tests);


end
