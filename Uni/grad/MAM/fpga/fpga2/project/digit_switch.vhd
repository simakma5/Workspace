LIBRARY ieee;
USE ieee.std_logic_1164.all;

ENTITY digit_switch IS
	PORT (
		CLK : IN STD_LOGIC;
		SW1 : IN STD_LOGIC;
		SW0 : IN STD_LOGIC;
		D0 : OUT STD_LOGIC;
		D1 : OUT STD_LOGIC;
		D2 : OUT STD_LOGIC;
		D3 : OUT STD_LOGIC;
		h : OUT STD_LOGIC
	);
END digit_switch;

ARCHITECTURE BEHAVIOR OF digit_switch IS
	TYPE digit IS (digit0,digit1,digit2,digit3);
	SIGNAL state : digit;
	SIGNAL next_state : digit;
	SIGNAL D : STD_LOGIC_VECTOR(3 downto 0);
BEGIN
	PROCESS (CLK,next_state)
	BEGIN
		IF(CLK='1' AND CLK'event) THEN
			state <= next_state;
		END IF;
	END PROCESS;

	PROCESS(state)
	BEGIN
		CASE state IS
			WHEN digit0 =>
			next_state <= digit1;
				D0 <= '1';
				D1 <= '0';
				D2 <= '0';
				D3 <= '0';
				D(0) <= '1';
				D(1) <= '0';
				D(2) <= '0';
				D(3) <= '0';
			WHEN digit1 =>
				next_state <= digit2;
				D0 <= '0';
				D1 <= '1';
				D2 <= '0';
				D3 <= '0';
				D(0) <= '0';
				D(1) <= '1';
				D(2) <= '0';
				D(3) <= '0';
			WHEN digit2 =>
				next_state <= digit3;
				D0 <= '0';
				D1 <= '0';
				D2 <= '1';
				D3 <= '0';
				D(0) <= '0';
				D(1) <= '0';
				D(2) <= '1';
				D(3) <= '0';
			WHEN digit3 =>
				next_state <= digit0;
				D0 <= '0';
				D1 <= '0';
				D2 <= '0';
				D3 <= '1';
				D(0) <= '0';
				D(1) <= '0';
				D(2) <= '0';
				D(3) <= '1';
			WHEN OTHERS => 
				D0 <= 'X';
				D1 <= 'X';
				D2 <= 'X';
				D3 <= 'X';
				D(0) <= 'X';
				D(1) <= 'X';
				D(2) <= 'X';
				D(3) <= 'X';
		END CASE;
	END PROCESS;
	 
	PROCESS(SW1, SW0, D)
	BEGIN
		IF (SW1 = '0') and (SW0 = '0') THEN
			IF D(1) = '1' THEN h <= '0';
			ELSE h <= '1'; 
			END IF;
		ELSIF (SW1 = '0') and (SW0 = '1') THEN
			IF D(2) = '1' THEN h <= '0';
			ELSE h <= '1';
			END IF;
		ELSIF (SW1 = '1') and (SW0 = '0') THEN
			h <= '1';
		ELSE h <= '1';
		END IF;
	END PROCESS;
END BEHAVIOR;
