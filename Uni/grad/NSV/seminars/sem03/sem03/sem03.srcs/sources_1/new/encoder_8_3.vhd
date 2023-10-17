library IEEE;
use IEEE.STD_LOGIC_1164.ALL;
use IEEE.STD_LOGIC_arith.ALL;
use IEEE.STD_LOGIC_unsigned.ALL;

entity encoder_8_3 is
    Port ( a : in STD_LOGIC_VECTOR (7 downto 0);
           b : out STD_LOGIC_VECTOR (2 downto 0));
end encoder_8_3;

architecture Behavioral of encoder_8_3 is

begin

process(a) begin
    b <= "000";
    for i in 0 to 7 loop
        if a(i) = '1' then
            b <= conv_std_logic_vector(i, 3);   --using a conversion function for convenience
        end if;
    end loop;
end process;


end Behavioral;
