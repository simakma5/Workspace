library IEEE;
use IEEE.STD_LOGIC_1164.ALL;

entity encoder_8_3_tb is
end encoder_8_3_tb;

architecture Behavioral of encoder_8_3_tb is

component encoder_8_3 is
    Port ( a : in STD_LOGIC_VECTOR (7 downto 0);
           b : out STD_LOGIC_VECTOR (2 downto 0));
end component;

signal a: std_logic_vector(7 downto 0) := "00000000";
signal b: std_logic_vector(2 downto 0);

begin

encoder_8_3_i: encoder_8_3
port map(a => a, b => b);

process begin
    a <= x"00";
    wait for 10ns;
    a <= x"01";
    wait for 10ns;
    a <= x"02";
    wait for 10ns;
    a <= x"04";
    wait for 10ns;
    a <= x"08";
    wait for 10ns;
    a <= x"0F";
    wait for 10ns;
    a <= x"AF";
    wait;
end process;

end Behavioral;
